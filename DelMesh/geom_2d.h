#ifndef GEOM_2D_H
#define GEOM_2D_H

#include "priorityqueue.h"

typedef Mesh::Point Point_2D;
class Geom_2D
{
public:
    Geom_2D(Mesh* mesh, Prop *samples);
    int get_sample_point(Mesh::EdgeHandle &ehandle);

    Mesh::Scalar get_interval(Mesh::HalfedgeHandle &heh, Mesh::Scalar &i1);

    void mesh_to_plane(Mesh::HalfedgeHandle heh, Mesh::Point &samp, Point_2D &p0, Point_2D &p1, Point_2D &dest);

    Mesh::Scalar distance2d(double *p1, Point_2D p2);

    Mesh::Scalar distance3d(Mesh::Point p1, Mesh::Point p2);

    void get_circumcircle(Point_2D &p0, Point_2D &p1, Point_2D &p2, double *pc, double &r);

    Mesh::Scalar square(Mesh::Scalar x);

    Point_2D get_2D_point(Mesh::HalfedgeHandle heh, Point_2D &p1, Point_2D &p2);

    vector<Point_2D> test_flattening(Mesh::EdgeHandle ehandle);

    void output_point(Mesh::EdgeHandle eh);


private:
    Mesh* mesh;
    Prop* samples;
};

#endif // GEOM_2D_H
