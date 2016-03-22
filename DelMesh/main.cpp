/*
 * delmesh.cpp
 *
 *  Created on: Mar 14, 2016
 *      Author: darryl
 */
#include "Delmesh.h"

int main(int argc, char** argv) {

    DelMesh mesh;
    //read the input mesh
    //mesh.test_pq();
    //mesh.test_face_colors();
    //mesh.test_2D_flattening();
    //mesh.test_samples();
    mesh.process_mesh();
    //mesh.test_face_colors();

    return 0;
}

