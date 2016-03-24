//============================================================================
// Name        : Delmesh.cpp
// Author      : darryl
// Version     :
// Copyright   : sure
// Description : Delaunay mesh generator
//============================================================================

#include <iostream>
#include "Delmesh.h"
#include <functional>
#include <math.h>
#include <OpenMesh/Core/Mesh/PolyMeshT.hh>

using namespace std;
#define INPUT "ateneav2.obj"
//#define INPUT "cube2.obj"
//#define OUTPUT "cube-out2.obj"
#define OUTPUT "ate-out.obj"
//0 = priority queue
//1 = queue
//2 = stack
#define DATA_STRUCTURE 0
#define score_type 2
#define TRUE 1
#define FALSE 0


DelMesh::DelMesh(){

    if (!OpenMesh::IO::read_mesh(mesh, INPUT)){
        cout << "read error"<<endl;
    }


    // Add the samples property, which is a vector of potential vertex
    // sites along an edge. See associated paper.
    mesh.add_property(samples);

    //add indicator variable for non-Delaunay edges
    mesh.add_property(is_NDE);

    //add indicator for a flippable edge
    mesh.add_property(is_flippable);

    //requesting face colors, although currently not working
    mesh.request_face_colors();
    mesh.request_face_status();
    mesh.request_edge_status();
    mesh.request_vertex_status();

    //make the constants we need
    make_constants();

    cout<<"pv: "<<pv<<endl;
    cout<<"pe: "<<pe<<endl;

    //non-Delaunay edges data structure

    q = new my_p_queue(&mesh, DATA_STRUCTURE);
    g2d = new Geom_2D(&mesh, &samples, &is_flippable, &is_NDE, score_type);


}

DelMesh::~DelMesh(){
    //write the mesh to file
    cout<<"writing output.obj"<<endl;
    if (!OpenMesh::IO::write_mesh(mesh, OUTPUT, writeOptions)){
        cout << "write error"<< endl;
    }
    delete q;
    delete g2d;
}

void DelMesh::process_mesh(){
    // find all the NDE's
    find_nd_edges();

    cout<<"done finding nde's"<<endl;

    make_Delaunay_mesh();

    sanity_check();

    while(q->size()>0){
        make_Delaunay_mesh();
        sanity_check();
    }

}

void DelMesh::sanity_check(){
    Mesh::EdgeIter eBegin, eEnd, eIt;
    eBegin = mesh.edges_begin();
    eEnd = mesh.edges_end();

    int count = 0;

    for (eIt = eBegin; eIt != eEnd; eIt++){
        //add sample points to the edge
        if (!mesh.is_boundary(*eIt)){
            if (is_nd_edge(*eIt, false)){
                if (mesh.property(is_flippable, *eIt)==FALSE){
                    if (mesh.property(samples, *eIt).size()==0){
                        cout << "bad edge, length "<<mesh.calc_edge_length(*eIt)<<endl;
                        count ++;
                    }else{
                        q->push(*eIt);
                    }

                }else{
                    cout<<"flipping"<<endl;
                    mesh.flip(*eIt);
                }
            }
        }
    }
    cout<<"Non-Delaunay edges: "<<count<<endl;
}

/**
     * @brief make_Delaunay_mesh
     * All NDE's have been found, now we process them and create the Delaunay mesh.
     */
void DelMesh::make_Delaunay_mesh(){

    Mesh::EdgeHandle eh;
    vector<Mesh::Point> samps;
    int index;
    Mesh::Point point;
    Mesh::Point p1, p2;
    Mesh::HalfedgeHandle heh1, heh2;
    Mesh::VertexHandle to, from, vh1, vh2, mid, v;
    Mesh::FaceHandle fh1, fh2;
    vector<Mesh::VertexHandle> vec;

    int dummy_count = 0;
    int count = 0;
    while(!q->empty()){
        //cout<<"counting iterations..."<<endl;
        count ++;
        //if (count >1) return;
        cout<<"edges: "<<mesh.n_edges()<<" count: "<<count<<" Stack size: "<<q->size()<<endl;
        eh = q->top();
        q->pop();

        //find the right sample point
        index = g2d->get_sample_point(eh);

        if (index == -1){
            //well fuck. There is no provision for this, because theoretically it should
            //never happen. But it does.
            continue;
            //we pretend it doesn't
        }

        //assuming that this invokes the copy constructor, since
        //we will be deleting the edge
        samps = mesh.property(samples, eh);

        point = samps.at(index);

        //cout<<"samps index point: "<<samps[index]<<endl;

        heh1 = mesh.halfedge_handle(eh, 0);
        heh2 = mesh.halfedge_handle(eh, 1);

        to = mesh.to_vertex_handle(heh1);
        from = mesh.from_vertex_handle(heh1);

        vh1 = mesh.opposite_vh(heh1);
        vh2 = mesh.opposite_vh(heh2);


        //grab the two faces we wish to delete
        fh1 = mesh.face_handle(heh1);
        fh2 = mesh.face_handle(heh2);

        mesh.delete_face(fh1, false);
        mesh.delete_face(fh2, false);
        mesh.garbage_collection();

        /*
         * Add the new vertex and the new faces to the mesh, while simultaneously
         * marking flippable edges, splitting up sample points, and checking and
         * adding new edges
         */

        //TODO: we should check the NDE flag before we add it a second time.
        // Flip flippable edges.
        mid = mesh.add_vertex(point);

        /* From here we add 4 faces. However, we need access to the edges, and the
         * only way to get them is using an iterator, which does not start at
         * any particular edge. So we have to test each edge using, in this case,
         * vertex handles.
         */

        Mesh::EdgeHandle eh;

        /************* Face 1 *******************/
        vec.clear();
        vec.push_back(from);
        vec.push_back(mid);
        vec.push_back(vh1);
        Mesh::FaceHandle fh = mesh.add_face(vec);

        //add in the samples in the split edge
        Mesh::FaceHalfedgeIter fhIt = mesh.fh_iter(fh);


        //iterate over the edges, figure out what we got and how to handle it
        for (int i = 0; i < 3; i ++){
            eh = mesh.edge_handle(*fhIt);
            v = mesh.from_vertex_handle(*fhIt);

            // edge (from, mid)
            if (v == from){
                cout<< "from,mid, vh1: v ==from "<<i<<endl;

                //we want to put the samples on, but we have to
                //put the correct ones
                if (is_on_edge(samps[0], eh)){
                    for (int i = 0; i < index; i++){
                        mesh.property(samples, eh).push_back(samps[i]);
                    }
                }else if (is_on_edge(samps[samps.size()-1], eh)){
                    for (int i = index +1; i < samps.size(); i++){
                        mesh.property(samples, eh).push_back(samps[i]);
                    }
                }else{
                    //put a dummy node for testing
                    cout<<"**********DUMMY NODE 1**************"<<endl;
                    mesh.property(samples, eh).push_back(Mesh::Point(1000,0,0));
                    dummy_count ++;
                }
                int c = mesh.property(samples, eh).size();
                cout<<c<<" sample edge 1******************************************"<<endl;
                // edge (mid, vh1)
            }else if (v == mid){
                cout<< "from,mid, vh1: v ==mid "<<i<<endl;

                mesh.property(is_flippable, eh) = TRUE;
                // edge (vh1, from)
            }else if (v == vh1){
                cout<< "from,mid, vh1: v ==vh1 "<<i<<endl;
                if (mesh.property(is_NDE, eh)==1){
                    if (mesh.property(is_flippable, eh)==TRUE){
                        test_flip(eh);
                        //mesh.flip(eh);

                    }
                }else if (is_nd_edge(eh, false)){
                    if (mesh.property(is_flippable, eh)==TRUE){
                        test_flip(eh);
                        //mesh.flip(eh);
                    }else{
                        mesh.property(is_NDE, eh) = TRUE;
                        //cout<<"pushing vh1, from"<<endl;
                        //g2d->output_point(eh);
                        q->push(eh);
                    }
                }

            }else{
                cout<<"damn"<<endl;
            }
            fhIt++;
        }

        /************* Face 2 *******************/
        vec.clear();
        vec.push_back(from);
        vec.push_back(vh2);
        vec.push_back(mid);
//        cout<<"from: "<<mesh.point(from)<<endl;
//        cout<<"vh2: "<<mesh.point(vh2)<<endl;
//        cout<<"mid: "<<mesh.point(mid)<<endl;

        fh = mesh.add_face(vec);
        fhIt = mesh.fh_iter(fh);

        //iterate over the edges, figure out what we got and how to handle it
        //we have already added samples above, so we skip that step
        for (int i = 0; i < 3; i ++){
            eh = mesh.edge_handle(*fhIt);
            v = mesh.from_vertex_handle(*fhIt);
            if (v==mid){
                int c = mesh.property(samples, eh).size();
                if (c ==0){cout<<" 0 sample edge 2****************************************** "<<i<<endl;}
                else{cout<<c<<" sample edge 2****************************************** "<<i<<endl;}
            }

            // edge (mid, from)
            if (v == vh2){
                cout<< "v ==vh2 "<<i<<endl;
                mesh.property(is_flippable, eh) = TRUE;
            // edge (from, vh2)
            }else if (v == from){
                cout<< "v ==from "<<i<<endl;
                if (mesh.property(is_NDE, eh)==TRUE){
                    if (mesh.property(is_flippable, eh)==TRUE){
                        test_flip(eh);
                        //mesh.flip(eh);
                    }
                }else if (is_nd_edge(eh, false)){
                    if (mesh.property(is_flippable, eh)==TRUE){
                        test_flip(eh);
                        //mesh.flip(eh);
                    }else{
                        mesh.property(is_NDE, eh) = TRUE;
                        //cout<<"pushing from, vh2"<<endl;
                        //g2d->output_point(eh);
                        q->push(eh);
                    }
                }

            }
            fhIt++;
        }

        /************* Face 3 *******************/
        vec.clear();
        vec.push_back(to);
        vec.push_back(vh1);
        vec.push_back(mid);

        fh = mesh.add_face(vec);
        fhIt = mesh.fh_iter(fh);

        //iterate over the edges, figure out what we got and how to handle it
        for (int i = 0; i < 3; i ++){
            eh = mesh.edge_handle(*fhIt);
            v = mesh.from_vertex_handle(*fhIt);

            // edge (mid, to)
            if (v == mid){
                //we want to put the samples on, but we have to
                //put the correct ones
                if (is_on_edge(samps[0], eh)){
                    for (int i = 0; i < index; i++){
                        mesh.property(samples, eh).push_back(samps[i]);
                    }
                }else if (is_on_edge(samps[samps.size()-1], eh)){
                    for (int i = index +1; i < samps.size(); i++){
                        mesh.property(samples, eh).push_back(samps[i]);
                    }
                }else{
                    //put a dummy node for testing
                    cout<<"**********DUMMY NODE 2**************"<<endl;
                    cout<<"samples "<<samps.size()<<" index: "<<index<<endl;
                    mesh.property(samples, eh).push_back(Mesh::Point(1000,0,0));
                    dummy_count ++;
                }
                int c = mesh.property(samples, eh).size();
                if (c ==0){cout<<" 0 sample edge 3 ******************************************"<<endl;}
            // edge (vh1, mid) has been handled
            // edge (to, vh1)
            }else if (v == to){
                if (mesh.property(is_NDE, eh)==TRUE){
                    if (mesh.property(is_flippable, eh)==TRUE){
                        test_flip(eh);
                        //mesh.flip(eh);
                    }
                }else if (is_nd_edge(eh, false)){
                    if (mesh.property(is_flippable, eh)==TRUE){
                        test_flip(eh);
                        //mesh.flip(eh);
                    }else{
                        mesh.property(is_NDE, eh) = TRUE;
                        //cout<<"pushing to, vh1"<<endl;
                        //g2d->output_point(eh);
                        q->push(eh);
                    }
                }

            }
            fhIt++;
        }

        /************* Face 4 *******************/
        vec.clear();
        vec.push_back(mid);
        vec.push_back(vh2);
        vec.push_back(to);
        fh = mesh.add_face(vec);
        fhIt = mesh.fh_iter(fh);

        //iterate over the edges, figure out what we got and how to handle it
        for (int i = 0; i < 3; i ++){
            eh = mesh.edge_handle(*fhIt);
            v = mesh.from_vertex_handle(*fhIt);

            if (v==to){
                int c = mesh.property(samples, eh).size();
                if (c ==0){cout<<" 0 sample edge ******************************************"<<endl;}
            }

            // edge (to, mid) has been handled
            // edge (mid, vh2) has been handled
            // edge (vh2, to)
            if (v == vh2){
                if (mesh.property(is_NDE, eh)==TRUE){
                    if (mesh.property(is_flippable, eh)==TRUE){
                        test_flip(eh);
                        //mesh.flip(eh);
                    }
                }else if (is_nd_edge(eh, false)){
                    if (mesh.property(is_flippable, eh)==TRUE){
                        test_flip(eh);
                        //mesh.flip(eh);
                    }else{
                        mesh.property(is_NDE, eh) = TRUE;
                        //cout<<"pushing vh2, to"<<endl;
                        //g2d->output_point(eh);
                        q->push(eh);
                    }
                }

            }
            fhIt++;
        }

        //cout<<"stack size: "<<q->size()<<endl;

        //I think that is everything. Looks like bug heaven.

        /*
         * find the opposite vertices, get the handles.
         * mark the two faces for deletion
         * call the garbage collector (probably)
         * add the sample vertex to the mesh
         * add the four new faces to the mesh
         *
         * We then have to check all 6 edges, two we have added and - well
         * We will grab face_handles for all the faces, iterate over the edges
         * to see if they are NDE, and add them to the q as necessary.
         *
         */
    }
    cout<<"dummy nodes: "<<dummy_count<<endl;

}

void DelMesh::test_flip(Mesh::EdgeHandle eh){
//    Mesh::VertexHandle to, from, t, f;
//    Mesh::HalfedgeHandle heh = mesh.halfedge_handle(eh, 0);

//    to = mesh.to_vertex_handle(heh);
//    from = mesh.from_vertex_handle(heh);
    mesh.flip(eh);

//    Mesh::HalfedgeHandle heh2 = mesh.halfedge_handle(eh, 0);
//    t = mesh.to_vertex_handle(heh2);
//    f = mesh.from_vertex_handle(heh2);


}

/**
 * We cheat a bit here. We find the midpoint of the line. Then see if the
 * vertex is within half the length of the midpoint.
 * @brief DelMesh::is_on_line
 * @param v
 * @param eh
 * @return
 */
bool DelMesh::is_on_edge(Mesh::Point &v, Mesh::EdgeHandle &eh)
{
    Mesh::Scalar length = mesh.calc_edge_length(eh)/2;
    Mesh::HalfedgeHandle heh = mesh.halfedge_handle(eh, 0);
    Mesh::Point to, from;
    to = mesh.point(mesh.to_vertex_handle(heh));
    from = mesh.point(mesh.from_vertex_handle(heh));
    if ((equals(v,to))||(equals(v,from))){
        return false;
    }

    to -= from;
    to.normalize();
    to *= length;
    to += from;
    bool answer = (g2d->distance3d(to, v)<= length);
//    if (answer){
//        cout<<"on edge"<<endl;
//    }else{
//        cout<<"NOT on edge"<<endl;
//    }
    return answer;
}

bool DelMesh::equals(Mesh::Point p1, Mesh::Point p2){
    bool answer = ((p1[0]==p2[0])&&(p1[1]==p2[1])&&(p1[2]==p2[2]));
//    if (answer){
//        cout<<"we have a sample point equal to the end point"<<endl;
//    }
    return answer;
}



void DelMesh::make_constants(){
    /*
         * Find the shortest edge, and the minimum angle
         */
    Mesh::VertexIter vIt, vBegin, vEnd;
    Mesh::VertexOHalfedgeIter veIt, veBegin, veEnd;;

    Mesh::Scalar min_angle = 7;
    Mesh::Scalar min_edge = 10000000;

    vBegin = mesh.vertices_begin();
    vEnd = mesh.vertices_end();
    int count = 0, c2 = 0;
    Mesh::Scalar length = 0, angle,maxl = 0;
    for (vIt = vBegin; vIt != vEnd; ++vIt){
        veBegin = mesh.voh_begin(*vIt);
        veEnd = mesh.voh_end(*vIt);

        for (veIt = veBegin; veIt != veEnd; ++veIt){
            c2++;
            length = mesh.calc_edge_length(*veIt);
            if (length > maxl){maxl = length;}
            if (length<min_edge && length != 0){
                min_edge = length;
            }
            angle = mesh.calc_sector_angle(*veIt);
            if (angle < 0){angle = - angle;}
            if (angle < min_angle&&angle != 0){
                min_angle = angle;
            }
            if (angle ==0){
                count++;
            }
        }
    }
    cout<<"zero angles: "<< count <<endl;
    cout<<"total angles: "<< c2 <<endl;
    cout<<"min angle: "<< min_angle <<endl;
    cout << "max length: "<< maxl<<endl;
    cout << "min length " << min_edge<<endl;


    pv = min_edge*sin(min_angle)/(0.5 + sin(min_angle));
    pe = min_edge/2;

    if (pe<pv){pv = pe;}

    pe = 2 * pv * sin(min_angle);

    double k = maxl/(min_edge*sin(min_angle)*sin(min_angle));

    cout <<"K: "<<k<<endl;

    double kk = maxl/ pe;

    cout <<"KK: "<<kk<<endl;
}


/*
     * We have been given a non-Delaunay edge, and now we will
     * add sample points to it.
     */
void DelMesh::make_sample_points(Mesh::EdgeHandle ehandle){
    /*
         * Get the two vertices associated with the edge, calculate and add samples
         * in a specified order. We always add the first sample point; after that
         * we check. For a definition of sample point, see the associated paper.
         */
    Mesh::HalfedgeHandle hedge = mesh.halfedge_handle(ehandle, 1);
    Mesh::Point to, from;
    to = mesh.point(mesh.to_vertex_handle(hedge));
    from = mesh.point(mesh.from_vertex_handle(hedge));

    double length = mesh.calc_edge_length(hedge);
    Mesh::Normal unit(mesh.calc_edge_vector(hedge));
    unit.normalize();

    double total = pv;
    Mesh::Point samp(unit);
    samp *= total;
    samp += from;
    mesh.property(samples, ehandle).push_back(samp);

    // We have our first sample in. If this is the shortest edge, it could
    // be the only sample (highly unlikely though).
    if (length <= 2 * pv){
        return;
    }

    // Otherwise keep adding samples
    double mark = length - (pv + pe);

    while (total < mark){
        total += pe;
        samp = (unit * total) + from;
        mesh.property(samples, ehandle).push_back(samp);
    }

    // We add one more sample at length *pv from the "to" vertex
    samp = (unit * (length - pv))+from;
    mesh.property(samples, ehandle).push_back(samp);

    // all samples have been added

}


/*
     * Run some tests to make sure the samples generated are accurate
     */
void DelMesh::test_samples(){
    Mesh::EdgeIter edge = mesh.edges_begin();
    make_sample_points(*edge);
    vector<Mesh::Point>* samps;
    cout<<"getting test samples"<<endl;
    samps = &(mesh.property(samples, *edge));
    cout<<"got test samples"<<endl;
    double length = mesh.calc_edge_length(*edge);
    int num_samples = 2;
    length -= 2*pv;

    length /= pe;
    num_samples += floor(length);
    cout << "expected samples: "<<num_samples<<endl;
    cout << "samples generated: "<<samps->size()<<endl;

    //now we want to check that they are all on the same line.
    //They should all output the same unit vector as the edge itself
    Mesh::Normal unit(mesh.calc_edge_vector(*edge));
    unit.normalize();
    cout<<endl<<"Normal vector: "<<unit<<endl;

    Mesh::HalfedgeHandle hedge = mesh.halfedge_handle(*edge, 1);

    Mesh::Point from = mesh.point(mesh.from_vertex_handle(hedge));
    Mesh::Point point(0,0,0);
    vector<Mesh::Point>::iterator it, it2;
    for (it = samps->begin(); it < samps->end(); it++){
        point = *it - from;
        point.normalize();
        cout<<endl<<"Normal vector: "<<unit<<endl;
        cout<<"Sample vector: "<<point<<endl;
        cout<<"Before normal: "<<*it<<endl;
    }

    edge ++; edge ++; edge ++;
    make_sample_points(*edge);
    vector<Mesh::Point>* sam;
    sam = &(mesh.property(samples, *edge));
    it = samps->begin(); it2 = sam->begin();
    for (int i = 0; i < 10; i++){
        cout << "Point 1: "<<*it++<<endl;
        cout << "Point 2: "<<*it2++<<endl;

    }


}



/**
     * We check the angles opposite the edge to see if they
     * sum to > pi.
     */
bool DelMesh::is_nd_edge(Mesh::EdgeHandle edge, bool output = false){

    Mesh::HalfedgeHandle hedge = mesh.halfedge_handle(edge, 1);
    Mesh::HalfedgeLoopIter hIt = mesh.hl_begin(hedge);

    //we want the angle opposite this half-edge, which is the
    //sector angle of the next half-edge
    hIt++;
    Mesh::Scalar angle1 = mesh.calc_sector_angle(*hIt);

    // we need absolute values
    if (angle1 <0){
        angle1 = -angle1;
    }

    //repeat for the opposite half-edge to get the opposite angle
    Mesh::HalfedgeHandle he_opposite_handle = mesh.opposite_halfedge_handle(hedge);
    hIt = mesh.hl_begin(he_opposite_handle);
    hIt++;
    Mesh::Scalar angle2 = mesh.calc_sector_angle(*hIt);
    if (angle2 <0){
        angle2 = -angle2;
    }
    if (output&&(angle1+angle2)>M_PI){
        cout<<"angle: "<<angle1+angle2<<endl;
    }
    return ((angle1+angle2)>M_PI);

}
/*
     * Find all non-Delaunay edges and put them in a priority queue for
     * max length
     */
void DelMesh::find_nd_edges(){
    Mesh::EdgeIter eBegin, eEnd, eIt;
    eBegin = mesh.edges_begin();
    eEnd = mesh.edges_end();


    // I will have to come back to this. The danger is that this edge might
    // be on a boundary. I think there is a test for that.
    //TODO boundary test
    for (eIt = eBegin; eIt != eEnd; eIt++){
        //add sample points to the edge
        if (!mesh.is_boundary(*eIt)){
            make_sample_points(*eIt);
            if (is_nd_edge(*eIt)){
                //bullshit doesn't work
                //mesh.set_color(*eIt, Mesh::Color(0.0,0.0,1.0,1.0));

                //indicate that the edge is non-Delaunay
                mesh.property(is_NDE, *eIt) = 1;
                //add the edge to our data structure of current NDE's
                q->push(*eIt);
            }
        }
    }
    cout<<"# NDE's: "<<q->size()<<endl;

}

Mesh* DelMesh::getMesh(){
    return &mesh;
}

Prop* DelMesh::getSamples(){
    return &samples;
}


/*
     * Test the priority queue. Presumably the NDE's have been added. Now
     * pop them off and output the lengths
     */
void DelMesh::test_pq(){
    Mesh::EdgeIter eBegin, eEnd, eIt;
    eBegin = mesh.edges_begin();
    eEnd = mesh.edges_end();

    /*
         * test it on 100 edges
         */

    Mesh::Scalar length = 0;
    for (int i = 0; i < 100; i++, eIt++){
        length = mesh.calc_edge_length(q->top());
        q->pop();
        cout << "length "<<i<<": "<<length<<endl;
    }

}

/*
     * Testing the face colors. Spoiler alert, they don't work. That is when
     * I try to load the file into a viewer after it throws and error.
     */
void DelMesh::test_face_colors(){
    if (!OpenMesh::IO::read_mesh(mesh, "ateneav.obj")){
        cout << "read error"<<endl;
    }
    if (!OpenMesh::IO::write_mesh(mesh, "ateneav2.obj", writeOptions)){
        cout << "write error"<< endl;
    }

}


int DelMesh::test_2D_flattening(){
    Mesh mesh;
    // generate vertices
    Mesh::VertexHandle vhandle[8];
    vhandle[0] = mesh.add_vertex(Mesh::Point(-1, -1,  1));
    vhandle[1] = mesh.add_vertex(Mesh::Point( 1, -1,  1));
    vhandle[2] = mesh.add_vertex(Mesh::Point( 1,  1,  1));
    vhandle[3] = mesh.add_vertex(Mesh::Point(-1,  1,  1));
    vhandle[4] = mesh.add_vertex(Mesh::Point(-1, -1, -1));
    vhandle[5] = mesh.add_vertex(Mesh::Point( 1, -1, -1));
    vhandle[6] = mesh.add_vertex(Mesh::Point( 1,  1, -1));
    vhandle[7] = mesh.add_vertex(Mesh::Point(-1,  1, -1));
    // generate (triangle) faces
    std::vector<Mesh::VertexHandle>  face_vhandles;
    face_vhandles.clear();
    face_vhandles.push_back(vhandle[0]);
    face_vhandles.push_back(vhandle[1]);
    face_vhandles.push_back(vhandle[2]);
    mesh.add_face(face_vhandles);


    face_vhandles.clear();
    face_vhandles.push_back(vhandle[0]);
    face_vhandles.push_back(vhandle[2]);
    face_vhandles.push_back(vhandle[3]);
    mesh.add_face(face_vhandles);

    //*********************/
    face_vhandles.clear();
    face_vhandles.push_back(vhandle[7]);
    face_vhandles.push_back(vhandle[6]);
    face_vhandles.push_back(vhandle[5]);
    mesh.add_face(face_vhandles);

    face_vhandles.clear();
    face_vhandles.push_back(vhandle[4]);
    face_vhandles.push_back(vhandle[7]);
    face_vhandles.push_back(vhandle[5]);
    mesh.add_face(face_vhandles);

    /**************/
    face_vhandles.clear();
    face_vhandles.push_back(vhandle[1]);
    face_vhandles.push_back(vhandle[0]);
    face_vhandles.push_back(vhandle[4]);
    mesh.add_face(face_vhandles);

    face_vhandles.clear();
    face_vhandles.push_back(vhandle[1]);
    face_vhandles.push_back(vhandle[4]);
    face_vhandles.push_back(vhandle[5]);
    mesh.add_face(face_vhandles);

    /*******************/
    face_vhandles.clear();
    face_vhandles.push_back(vhandle[2]);
    face_vhandles.push_back(vhandle[5]);
    face_vhandles.push_back(vhandle[6]);
    mesh.add_face(face_vhandles);

    face_vhandles.clear();
    face_vhandles.push_back(vhandle[1]);
    face_vhandles.push_back(vhandle[5]);
    face_vhandles.push_back(vhandle[2]);
    mesh.add_face(face_vhandles);

    /****************/
    face_vhandles.clear();
    face_vhandles.push_back(vhandle[2]);
    face_vhandles.push_back(vhandle[6]);
    face_vhandles.push_back(vhandle[3]);
    mesh.add_face(face_vhandles);

    face_vhandles.clear();
    face_vhandles.push_back(vhandle[3]);
    face_vhandles.push_back(vhandle[6]);
    face_vhandles.push_back(vhandle[7]);
    mesh.add_face(face_vhandles);

    /**************/
    face_vhandles.clear();
    face_vhandles.push_back(vhandle[0]);
    face_vhandles.push_back(vhandle[3]);
    face_vhandles.push_back(vhandle[7]);
    mesh.add_face(face_vhandles);

    face_vhandles.clear();
    face_vhandles.push_back(vhandle[0]);
    face_vhandles.push_back(vhandle[7]);
    face_vhandles.push_back(vhandle[4]);
    mesh.add_face(face_vhandles);

    // write mesh to cube.obj

    try
    {
        if ( !OpenMesh::IO::write_mesh(mesh, "cube.obj") )
        {
            std::cerr << "Cannot write mesh to file 'cube.obj'" << std::endl;
        }
    }
    catch( std::exception& x )
    {
        std::cerr << x.what() << std::endl;
    }


    Geom_2D g(&mesh, &samples, &is_flippable, &is_NDE, score_type);

    Mesh::EdgeIter eIt = mesh.edges_begin();

    Mesh::HalfedgeHandle hedge = mesh.halfedge_handle(*eIt, 1);
    Mesh::Point to, from;
    to = mesh.point(mesh.to_vertex_handle(hedge));
    from = mesh.point(mesh.from_vertex_handle(hedge));

    cout<<"point 1 x:"<<to[0]<<" y: "<<to[1]<<" z: "<<to[2]<<endl;
    cout<<"point 2 x:"<<from[0]<<" y: "<<from[1]<<" z: "<<from[2]<<endl;


    vector<Point_2D> p(g.test_flattening(*eIt));


    Mesh mesh2;
    // generate vertices
    for (int i = 0; i < 8; i++){
        vhandle[i] = mesh2.add_vertex( (p[i]));
        std::cout<<"point "<<i<<": "<<(p[i])<<endl;
    }

    face_vhandles.clear();
    face_vhandles.push_back(vhandle[0]);
    face_vhandles.push_back(vhandle[1]);
    face_vhandles.push_back(vhandle[2]);
    mesh2.add_face(face_vhandles);
    cout<<"complex?"<<endl;

    face_vhandles.clear();
    face_vhandles.push_back(vhandle[0]);
    face_vhandles.push_back(vhandle[5]);
    face_vhandles.push_back(vhandle[1]);
    mesh2.add_face(face_vhandles);
    cout<<"complex?"<<endl;

    /*********************/
    face_vhandles.clear();
    face_vhandles.push_back(vhandle[7]);
    face_vhandles.push_back(vhandle[1]);
    face_vhandles.push_back(vhandle[5]);
    mesh2.add_face(face_vhandles);
    cout<<"complex?"<<endl;

    face_vhandles.clear();
    face_vhandles.push_back(vhandle[0]);
    face_vhandles.push_back(vhandle[6]);
    face_vhandles.push_back(vhandle[5]);
    mesh2.add_face(face_vhandles);
    cout<<"complex?"<<endl;

    /**************/
    face_vhandles.clear();
    face_vhandles.push_back(vhandle[2]);
    face_vhandles.push_back(vhandle[1]);
    face_vhandles.push_back(vhandle[3]);
    mesh2.add_face(face_vhandles);
    cout<<"complex?"<<endl;

    face_vhandles.clear();
    face_vhandles.push_back(vhandle[0]);
    face_vhandles.push_back(vhandle[2]);
    face_vhandles.push_back(vhandle[4]);
    mesh2.add_face(face_vhandles);
    cout<<"complex?"<<endl;

    // write mesh to flat_cube.obj
    try
    {
        if ( !OpenMesh::IO::write_mesh(mesh2, "flat_cube.obj") )
        {
            std::cerr << "Cannot write mesh to file 'cube.obj'" << std::endl;
            return 1;
        }
    }
    catch( std::exception& x )
    {
        std::cerr << x.what() << std::endl;
        return 1;
    }
    return 0;
}








