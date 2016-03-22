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
//#define INPUT "ateneav.obj"
#define INPUT "cube2.obj"
#define OUTPUT "out-cube.obj"
//0 = priority queue
//1 = queue
//2 = stack
#define DATA_STRUCTURE 0


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
    g2d = new Geom_2D(&mesh, &samples);


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
    Mesh::VertexHandle to, from, vh1, vh2, mid;
    Mesh::FaceHandle fh1, fh2;
    vector<Mesh::VertexHandle> vec;

    int count = 0;
    while(!q->empty()){
        //cout<<"counting iterations..."<<endl;
        count ++;
        if (count >50) return;
        cout<<"edges: "<<mesh.n_edges()<<endl;
        eh = q->top();
        q->pop();

        //find the right sample point
        index = g2d->get_sample_point(eh);

        //assuming that this invokes the copy constructor, since
        //we will be deleting the edge
        samps = mesh.property(samples, eh);

        point = samps.at(index);

        cout<<"samps index point: "<<samps[index]<<endl;

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
        p1 = mesh.point(mid);
        cout <<"added point:" <<p1<<endl;

        vec.clear();
        vec.push_back(from);
        vec.push_back(mid);
        vec.push_back(vh1);
        Mesh::FaceHandle fh = mesh.add_face(vec);

        //add in the samples in the split edge
        Mesh::FECWWIter feIt = mesh.fe_ccwbegin(fh);
        for (int i = 0; i < index; i++){
            mesh.property(samples, *feIt).push_back(samps[i]);
        }

        //this is the new edge. It will always be Delaunay for this step.
        //There are no sample points, but it is flippable.
        feIt++;
        mesh.property(is_flippable, *feIt) = 1;

        //this edge was present before. It may be ND at this point.
        feIt++;
        if (mesh.property(is_NDE, *feIt)==1){
            if (mesh.property(is_flippable, *feIt)){
                mesh.flip(*feIt);
            }
        }else if (is_nd_edge(*feIt)){
            mesh.property(is_NDE, *feIt) = 1;
            cout<<"pushing1"<<endl;
            g2d->output_point(*feIt);
            q->push(*feIt);
        }

        vec.clear();
        vec.push_back(from);
        vec.push_back(vh2);
        vec.push_back(mid);
        cout<<"from: "<<mesh.point(from)<<endl;
        cout<<"vh2: "<<mesh.point(vh2)<<endl;
        cout<<"mid: "<<mesh.point(mid)<<endl;

        fh = mesh.add_face(vec);
        feIt = mesh.fe_ccwbegin(fh);
        //(to,vh2)
        if (mesh.property(is_NDE, *feIt)==1){
            if (mesh.property(is_flippable, *feIt)){
                mesh.flip(*feIt);
            }
        }else if (is_nd_edge(*feIt)){
            mesh.property(is_NDE, *feIt) = 1;
            cout<<"pushing2: "<<endl;
            g2d->output_point(*feIt);
            cout<<"edge length: "<<mesh.calc_edge_length(*feIt)<<endl;
            vector<Mesh::Point>* samps;
            samps = &(mesh.property(samples, *feIt));
            cout<<"sample size: "<<samps->size()<<endl;
            q->push(*feIt);
        }

        /*********************************testing****/
        feIt++;
        g2d->output_point(*feIt);
        feIt++;
        g2d->output_point(*feIt);

//        Mesh::Point p1 = mesh.point(mesh.to_vertex_handle(mesh.halfedge_handle(*feIt,0)));
//        Mesh::Point p2 = mesh.point(mesh.from_vertex_handle(mesh.halfedge_handle(*feIt,0)));

//        cout <<"point1: "<<p1<<endl;
//        cout <<"point2: "<<p2<<endl;


        vec.clear();
        vec.push_back(to);
        vec.push_back(vh1);
        vec.push_back(mid);
        fh = mesh.add_face(vec);
        feIt = mesh.fe_ccwbegin(fh);
        //(from, vh1)
        if (mesh.property(is_NDE, *feIt)==1){
            if (mesh.property(is_flippable, *feIt)){
                mesh.flip(*feIt);
            }
        }else if (is_nd_edge(*feIt)){
            mesh.property(is_NDE, *feIt) = 1;
            cout<<"pushing3"<<endl;
            g2d->output_point(*feIt);
            q->push(*feIt);
        }

        feIt++; //(vh1,mid)
        mesh.property(is_flippable, *feIt) = 1;
        //the other half of our split edge, (from, mid)
        feIt++;
        for (int i = index +1; i < samps.size(); i++){
            //cout<<"pushing back "<<i<<endl;
            mesh.property(samples, *feIt).push_back(samps[i]);
        }

        vec.clear();
        vec.push_back(mid);
        vec.push_back(vh2);
        vec.push_back(to);
        fh = mesh.add_face(vec);
        feIt = mesh.fe_ccwbegin(fh);
        mesh.property(is_flippable, *feIt) = 1; //(mid, vh2)

        feIt++; //(vh2, from)
        if (mesh.property(is_NDE, *feIt)==1){
            if (mesh.property(is_flippable, *feIt)){
                mesh.flip(*feIt);
            }
        }else if (is_nd_edge(*feIt)){
            mesh.property(is_NDE, *feIt) = 1;
            cout<<"pushing4"<<endl;

            q->push(*feIt);
        }

        cout<<"stack size: "<<q->size()<<endl;

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
bool DelMesh::is_nd_edge(Mesh::EdgeHandle edge){

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
        if (!mesh.is_boundary(*eIt)){
            if (is_nd_edge(*eIt)){
                //bullshit doesn't work
                //mesh.set_color(*eIt, Mesh::Color(0.0,0.0,1.0,1.0));

                //indicate that the edge is non-Delaunay
                mesh.property(is_NDE, *eIt) = 1;
                //add sample points to the edge
                make_sample_points(*eIt);
                //add the edge to our data structure of current NDE's
                q->push(*eIt);
            }
        }
    }

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


    Geom_2D g(&mesh, &samples);

    Mesh::EdgeIter eIt = mesh.edges_begin();

    Mesh::HalfedgeHandle hedge = mesh.halfedge_handle(*eIt, 1);
    Mesh::Point to, from;
    to = mesh.point(mesh.to_vertex_handle(hedge));
    from = mesh.point(mesh.from_vertex_handle(hedge));

    cout<<"point 1 x:"<<to[0]<<" y: "<<to[1]<<" z: "<<to[2]<<endl;
    cout<<"point 2 x:"<<from[0]<<" y: "<<from[1]<<" z: "<<from[2]<<endl;


    cout<<"got here"<<endl;
    vector<Point_2D> p(g.test_flattening(*eIt));

    cout<<"get here?"<<endl;

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








