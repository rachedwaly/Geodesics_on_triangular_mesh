#include <iostream>
#include <igl/readOFF.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <direct.h>
#include <igl/avg_edge_length.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/parula.h>
#include <igl/per_corner_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/principal_curvature.h>
#include <igl/barycenter.h>
#include <igl/grad.h>
#include <igl/jet.h>
#include <igl/internal_angles.h>
#include <igl/isolines_map.h>

#include "HalfedgeBuilder.cpp"



using namespace Eigen; // to use the classes provided by Eigen library
using namespace std;

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
    if (key == '1')
    {
        viewer.data().clear();
        
        
    }
    else if (key == '2')
    {
        viewer.data().clear();
        
    }
    return false;
}

MatrixXd areaMatrix(MatrixXd& V, MatrixXi& F) {
    MatrixXd res = MatrixXd::Zero(F.rows(), 1);
    Vector3d v1(0.0,0.0,0.0);
    Vector3d v2(0.0,0.0,0.0);
    Vector3d v3(0.0,0.0,0.0);
    Vector3d a, b;

    for (int i = 0; i < F.rows(); ++i) {
        v1 << V(F(i, 0), 0), V(F(i, 0), 1), V(F(i, 0), 2);
        v2 << V(F(i, 1), 0), V(F(i, 1), 1), V(F(i, 1), 2);
        v3 << V(F(i, 2), 0), V(F(i, 2), 1), V(F(i, 2), 2);
        a = v2 - v1;
        b = v3 - v1;
        res(i, 0) = (1.0 / 2) * (a.cross(b)).norm();
    }
    return res;
}
//this method will calculate the neighbour faces for a given vertex v
vector<int> vertexNeighbour(HalfedgeDS he, int v) {
    vector <int> n;
    int edge = he.getEdge(v);
    int oppEdge = he.getOpposite(edge);
    int prevEdge = he.getPrev(oppEdge);
    if (he.getFace(prevEdge) != -1) {
        n.push_back(he.getFace(prevEdge));
    }
    while (prevEdge != edge) {
        oppEdge = he.getOpposite(prevEdge);
        prevEdge = he.getPrev(oppEdge);
        if (he.getFace(prevEdge) != -1) {
            n.push_back(he.getFace(prevEdge));
        }
    }
    return n;
}



bool elementexists(vector<int> v, int j) {
    for (int i = 0; i < v.size(); i = i + 1) {
        if (v.at(i) == j) return true;
    }
    return false;
}

//this method will calculate the neighbour faces for a given face
vector<int> faceNeighbours(HalfedgeDS he, MatrixXd &V, MatrixXi& F,int f) {
    vector<int> res;
    for (int j = 0; j < 3; j = j + 1) {
            vector<int> a = vertexNeighbour(he, F(f, j));
            for (int i = 0; i < a.size(); i = i + 1) {
                if (!(elementexists(res, a.at(i)))) {
                    res.push_back(a.at(i));
                }
            }
        }
    return res;
}

int nextTriangleInColorFlow(HalfedgeDS he, MatrixXd& V, MatrixXi& F, int f) {
    return 0;
}


void buildNeighbours(std::map<int, vector<int>>& adj, HalfedgeDS he, MatrixXd& V) {
    for (int i = 0; i < V.rows(); i++) {
        adj.insert({ i, vertexNeighbour(he, i) });
    }
}




MatrixXd divMatrix1(std::map<int, vector<int>>& adj, MatrixXd& V,MatrixXi &F, MatrixXd& areaMatrix,MatrixXd & massMatrix, MatrixXd& grad) {
    MatrixXd div = MatrixXd::Zero(V.rows(), F.rows() * 3);

    for (int i = 0; i < V.rows(); i++) {
        
        for (const auto& value : adj.at(i)) {
                div(i, 3 * value) = -(areaMatrix(value, 0) / massMatrix(i, i)) * grad(3 * value, i);
                div(i, 3 * value + 1) = -(areaMatrix(value, 0) / massMatrix(i, i)) * grad(3 * value + 1, i);
                div(i, 3 * value + 2) = -(areaMatrix(value, 0) / massMatrix(i, i)) * grad(3 * value + 2, i);
        }
    }
    
    return div;
}

double distance(MatrixXd &V1, MatrixXd &V2) {
    return (V2 - V1).norm();
}

double cosinus_law(double a, double b, double c) {
    return acos((-pow(c, 2) + pow(a, 2) +pow(b, 2)) / (2.0 * a * b));
}

double cotan(double i) { return(1 / tan(i)); }

MatrixXd internal_angles(MatrixXd& V, MatrixXi& F) {
    MatrixXd angles(F.rows(), 3);
    double a=0, b=0, c=0;
    MatrixXd A(1,3), B(1,3), C(1,3);
    for (int i = 0; i < F.rows(); i = i + 1) {
        A = V.row(F(i, 0));
        B = V.row(F(i, 1));
        C = V.row(F(i, 2));

        a = distance(C, B);
        b = distance(A, C);
        c = distance(A, B);

        

        angles(i, 0) = cosinus_law(c, b, a);
        angles(i, 1) = cosinus_law(a, c, b);
        angles(i, 2) = cosinus_law(a, b, c);

    }

    return angles;
}

double cross(MatrixXd& e, MatrixXd& X) {
    Vector3d e1(e(0, 0), e(0, 1), e(0, 2));
    Vector3d X1(X(0, 0), X(0, 1), X(0, 2));
    return (e1.dot(X1));
}



MatrixXd divMatrix(std::map<int, vector<int>>& adj, MatrixXd& V, MatrixXi& F, MatrixXd& areaMatrix, MatrixXd& massMatrix, MatrixXd& gradU) {
    MatrixXd div = MatrixXd::Zero(V.rows(),1);
    MatrixXd angles = internal_angles(V, F);
   
    MatrixXd X(1, 3), e1(1, 3), e2(1, 3);
    vector<int> e,ind;
    
    for (int i = 0; i < V.rows(); i++) {
        
        for (const auto& value : adj.at(i)) {
            X = gradU.row(value);
            for (int j = 0; j < 3; j = j + 1) {
                if (F(value, j) != i) 
                    e.push_back(F(value, j));
                    ind.push_back(j);
            }
            e1 = V.row(e.at(0)) - V.row(i);
            e2 = V.row(e.at(1)) - V.row(i);
            div(i, 0) += 0.5*cotan(angles(value,ind.at(0)))*cross(e1,X)+0.5*cotan(angles(value,ind.at(1)))*cross(e2,X);
            
            e.clear();
            ind.clear();
        }
        }
    

    return div;
}
void set_colormap(igl::opengl::glfw::Viewer& viewer)
{
    const int num_intervals = 30;
    Eigen::MatrixXd CM(num_intervals, 3);
    // Colormap texture
    for (int i = 0; i < num_intervals; i++)
    {
        double t = double(num_intervals - i - 1) / double(num_intervals - 1);
        CM(i, 0) = std::max(std::min(2.0 * t - 0.0, 1.0), 0.0);
        CM(i, 1) = std::max(std::min(2.0 * t - 1.0, 1.0), 0.0);
        CM(i, 2) = std::max(std::min(6.0 * t - 5.0, 1.0), 0.0);
    }
    igl::isolines_map(Eigen::MatrixXd(CM), CM);
    viewer.data().set_colormap(CM);
}



void showVectorField(MatrixXd& V, MatrixXi& F, MatrixXd& GU, igl::opengl::glfw::Viewer& viewer) {
    // Draw a black segment in direction of gradient at face barycenters
    MatrixXd BC;
    igl::barycenter(V, F, BC);
    const RowVector3d red(1, 0, 0);
    const RowVector3d blue(0, 0, 1);
    viewer.data().add_edges(BC, BC + GU * 0.5, red);
    viewer.data().add_edges(BC + GU * 0.5, BC + GU * 0.5 + GU * 0.5, blue);
}

map<int, vector<int>> flowMapByFaces(MatrixXd& V, MatrixXi& F, MatrixXd& GU) {
    MatrixXd BC;
    map<int, vector<int>> mapOfFlowByfaces;
    igl::barycenter(V, F, BC);
    double c1, c2, c3;
    vector<int> a;
    double k = 0;
    Vector3d V1(0, 0, 0), V2(0, 0, 0), V3(0, 0, 0);
    for (int i = 0; i < F.rows(); i = i + 1) {
        V1 = V.row(F(i, 0)) - BC.row(i);
        V2 = V.row(F(i, 1)) - BC.row(i);
        V3 = V.row(F(i, 2)) - BC.row(i);
        /*
        if (V1.dot(GU.row(i)) > -0.1) {
            a.push_back(F(i, 0));
        }
        if (V2.dot(GU.row(i)) > -0.1) {
            a.push_back(F(i, 1));
        }
        if (V3.dot(GU.row(i)) > -0.1) {
            a.push_back(F(i, 2));
        }
        */
        c1 = acos(V1.dot(GU.row(i)) / (V1.norm() * GU.row(i).norm()));
        c2 = acos(V2.dot(GU.row(i)) / (V2.norm() * GU.row(i).norm()));
        c3 = acos(V3.dot(GU.row(i)) / (V3.norm() * GU.row(i).norm()));

        k=(max(c1, max(c2, c3)));
        
        if (k == c1) { 
            a.push_back(F(i, 2)); 
            a.push_back(F(i, 1));
            a.push_back(F(i, 0));
        }
        if (k == c2) { 
            a.push_back(F(i, 2)); 
            a.push_back(F(i, 0));
            a.push_back(F(i, 1));
        }
        if (k == c3) { 
            a.push_back(F(i, 0));
            a.push_back(F(i, 1));
            a.push_back(F(i, 2));
        }

        mapOfFlowByfaces.insert({i,vector<int>(a)});
        a.clear();
    }
    return mapOfFlowByfaces;
}

void showMaxMinCurvature(MatrixXd& V, MatrixXi& F, igl::opengl::glfw::Viewer& viewer) {
    MatrixXd PD1, PD2;
    VectorXd PV1, PV2;
    igl::principal_curvature(V, F, PD1, PD2, PV1, PV2);
    // mean curvature
    MatrixXd H = 0.5 * (PV1 + PV2);


    viewer.data().set_data(H);
    std::cout << PV1 << std::endl;
        // Average edge length for sizing
    const double avg = igl::avg_edge_length(V, F);

    // Draw a red segment parallel to the maximal curvature direction
    const RowVector3d red(0.8, 0.2, 0.2), blue(0.2, 0.2, 0.8);
    viewer.data().add_edges(V + PD1 * avg, V - PD1 * avg, red);

    // Draw a blue segment parallel to the minimal curvature direction
    viewer.data().add_edges(V + PD2 * avg, V - PD2 * avg, blue);
}


void color(map<int, vector<int>>& mapOfFlow,MatrixXd& V,int source,MatrixXd & c,MatrixXd & phi2,map<int,int> mapOfFlags) {
    
    map<int, vector<int>>::iterator it;
    
    for (it = mapOfFlow.begin(); it != mapOfFlow.end(); it++)
    {
        for (int i = 0; i < it->second.size() - 1; i = i + 1) {

            if (abs(phi2(it->second.at(i), 2) - phi2(it->second.at(2), 2)) > 0.15) {
                mapOfFlags.at(it->second.at(i)) = 1;
                c(it->second.at(i)) = 1;

            }

        
            
        }
    }
    

}

int lookForNext(int sourceID, MatrixXd V,MatrixXi F, MatrixXd minCurvD,HalfedgeDS he) {
    vector<int> neighbours = vertexNeighbour(he, sourceID);

    return 0;
}

int main()
{
    // Load a mesh in OFF format
    igl::opengl::glfw::Viewer viewer;
    Eigen::MatrixXd V, C;
    Eigen::MatrixXi F;

    Eigen::MatrixXd I,A;

    std::string filename = "C:\\Users\\Rached\\Documents\\Telecom\\igd\\Xinf574\\project\\data\\mctest2.off";
    igl::read_triangle_mesh(filename, V, F);


    viewer.data().set_mesh(V, F);
    
    HalfedgeBuilder* builder = new HalfedgeBuilder();  

    HalfedgeDS he = builder->createMesh(V.rows(), F);  
    
    std::map<int, vector<int>> adj;
    Vector3d qsdqs = V.row(0);


    std::cout << "calculating neighbours" << std::endl;
    buildNeighbours(adj, he, V);

    
    
    const double h = igl::avg_edge_length(V, F);
    double t = pow(h, 2);
    

    Eigen::MatrixXd U0=Eigen::MatrixXd::Zero(V.rows(), 1);




    Eigen::MatrixXd U1 = Eigen::MatrixXd::Zero(V.rows(), 1);
    


    
    U0(93, 0) = 1;
    
    
    

    
    
    
    Eigen::SparseMatrix<double> L, M;
    igl::cotmatrix(V, F, L);
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
    
    
    A = M - t * L;
    
    Eigen::FullPivLU<Eigen::MatrixXd> dec(A);
   
    U1 = dec.solve(U0);



    
    // Compute gradient operator:
    Eigen::SparseMatrix<double> G;
    igl::grad(V, F, G);

    
    
    
    Eigen::MatrixXd GU = Eigen::Map<const Eigen::MatrixXd>((G * U1).eval().data(), F.rows(), 3);


    
    MatrixXd div;



    MatrixXd area;
    std::cout << "calculating areas" << std::endl;
    area = areaMatrix(V, F);
    

    GU = -GU;
    GU.rowwise().normalize();


    
    cout << "Calculating geodesic paths" << endl;
    auto start = std::chrono::high_resolution_clock::now();
    
    div = divMatrix(adj, V, F, area, MatrixXd(M), GU);




    MatrixXd phi=MatrixXd::Zero(V.rows(),1);



    Eigen::FullPivLU<Eigen::MatrixXd> dec1(L);
    

    phi = dec1.solve(div);



    // for measuring time performances
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "calculation time: " << elapsed.count() << " s\n";
    
      
    //showVectorField(V, F, GU, viewer);
   
    
    // Initialize white
    C = Eigen::MatrixXd::Constant(F.rows(), 3, 1);
    //viewer.data().set_colors(C);
    map<int, vector<int>> mapOfFlow = flowMapByFaces(V, F, GU);

    MatrixXd colors(V.rows(), 1);

    viewer.callback_mouse_down = [&V, &F, &C,&he,&colors,&mapOfFlow](igl::opengl::glfw::Viewer& viewer, int, int)->bool
    {
        int fid;
        Eigen::Vector3f bc;
        // Cast a ray in the view direction starting from the mouse position
        double x = viewer.current_mouse_x;
        
       
        double y = viewer.core().viewport(3) - viewer.current_mouse_y;

        if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core().view,
            viewer.core().proj, viewer.core().viewport, V, F, fid, bc))
        {
            // paint hit red
            C.row(fid) << 1, 0, 0;
            
            /*
            std::cout << "index of face is " << fid << std::endl;
            std::cout << "vertices are " << F.row(fid) << std::endl;
            std::cout << phi.row(F(fid, 0)) << std::endl;
            std::cout << phi.row(F(fid, 1)) << std::endl;
            std::cout << phi.row(F(fid, 2)) << std::endl;
            

            for (const auto& value : vertexNeighbour(he, F(fid, 0))) {
                C.row(value) << 1, 0, 0;
            }
            for (const auto& value : vertexNeighbour(he, F(fid, 1))) {
                C.row(value) << 0, 1, 0;
            }
            for (const auto& value : vertexNeighbour(he, F(fid, 2))) {
                C.row(value) << 0, 0, 1;
            }
            */
            //viewer.data().set_colors(C);
            int sourceID = F(fid, 0);
            for (int i = 0; i < 10; i++) {
                //sourceID = lookForNext(sourceID, V,F, PD2.row(sourceID));
            }


            return true;
        }
        return false;
    };
        
    viewer.data().set_data(colors); 
    viewer.callback_key_down = &key_down;
    viewer.data().show_lines = true;
    viewer.launch();
    
}