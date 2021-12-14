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
#include <queue>

#include "HalfedgeBuilder.cpp"



using namespace Eigen; // to use the classes provided by Eigen library
using namespace std;

struct Point {
	int vertexIndex;
	double distanceToSource;
};

struct pointComp {
	bool operator() (Point a, Point b) {
		return a.distanceToSource > b.distanceToSource;
	}
};

MatrixXd phi2;

MatrixXd phi;

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
	if (key == '1')
	{
		viewer.data().set_data(phi);
		std::cout << "!!!!!!!!!!!!!!!! heat method !!!!!!!!!!!!!!!!" << std::endl;
	}
	else if (key == '2')
	{
		viewer.data().set_data(phi2);
		std::cout << "!!!!!!!!!!!!!!!! disjktra method !!!!!!!!!!!!!!!!" << std::endl;
	}
	return false;
}

MatrixXd areaMatrix(MatrixXd& V, MatrixXi& F) {
	MatrixXd res = MatrixXd::Zero(F.rows(), 1);
	Vector3d v1(0.0, 0.0, 0.0);
	Vector3d v2(0.0, 0.0, 0.0);
	Vector3d v3(0.0, 0.0, 0.0);
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



MatrixXd dijkstra(MatrixXd& V, MatrixXi& F, VectorXi sources) {

	int n = V.rows();
	MatrixXd phi = MatrixXd::Zero(n + 1, 1);
	MatrixXd dist = MatrixXd::Zero(n + 1, 1);
	MatrixXi visited = MatrixXi::Zero(n + 1, 1);


	dist(n, 0) = 0;

	priority_queue<Point, vector<Point>, pointComp> pq;

	Point startPoint({ n,0 });
	pq.push(startPoint);



	vector<int>* adjList = new  vector<int>[n + 2];
	for (int i = 0; i < n; i++) {
		dist(i, 0) = (double)INFINITY;
		if (sources(i, 0) == 1) {
			adjList[n].push_back(i);
		}
	}
	for (int i = 0; i < F.rows(); i++) {
		for (int j = 0; j < F.cols(); j++) {
			adjList[F(i, j)].push_back(F(i, (j + 1) % 3));
			adjList[F(i, j)].push_back(F(i, (j + 2) % 3));
		}
	}

	while (!pq.empty()) {
		int index;
		double currentDist;

		Point currentVertex = pq.top();
		pq.pop();
		index = currentVertex.vertexIndex;
		currentDist = currentVertex.distanceToSource;
		//cout << " # " << index << " # " << currentDist << endl;
		if (visited(index, 0) == 0) {
			phi(index, 0) = currentDist;
			visited(index, 0) = 1;

			for (int i = 0; i < adjList[index].size(); i++) {
				int neighborIndex = adjList[index][i];
				if (visited(neighborIndex, 0) == 0) {
					double edgeNorm = (index == n) ? 0 : (V.row(index) - V.row(neighborIndex)).norm();
					dist(neighborIndex, 0) = min(dist(neighborIndex, 0), currentDist + edgeNorm);
					Point a({ neighborIndex,  dist(neighborIndex, 0) });
					pq.push(a);
					//cout << neighborIndex << " : " << dist(neighborIndex, 0) << endl;
				}
			}
		}
	}

	delete[] adjList;
	return phi.topRows(n);
}

bool elementexists(vector<int> v, int j) {
	for (int i = 0; i < v.size(); i = i + 1) {
		if (v.at(i) == j) return true;
	}
	return false;
}

//this method will calculate the neighbour faces for a given face
vector<int> faceNeighbours(HalfedgeDS he, MatrixXd& V, MatrixXi& F, int f) {
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




double distance(MatrixXd& V1, MatrixXd& V2) {
	return (V2 - V1).norm();
}

double cosinus_law(double a, double b, double c) {
	return acos((-pow(c, 2) + pow(a, 2) + pow(b, 2)) / (2.0 * a * b));
}

double cotan(double i) { return(1 / tan(i)); }

MatrixXd internal_angles(MatrixXd& V, MatrixXi& F) {
	MatrixXd angles(F.rows(), 3);
	double a = 0, b = 0, c = 0;
	MatrixXd A(1, 3), B(1, 3), C(1, 3);
	for (int i = 0; i < F.rows(); i = i + 1) {
		A = V.row(F(i, 0));
		B = V.row(F(i, 1));
		C = V.row(F(i, 2));

		a = distance(C, B);
		b = distance(A, C);
		c = distance(A, B);

		angles(i, 0) = cosinus_law(c, b, a);
		angles(i, 1) = cosinus_law(a, c, b);
		angles(i, 2) = cosinus_law(b, a, c);

	}

	return angles;
}

double dot(const MatrixXd& e, const MatrixXd& X) {
	Vector3d e1(e(0, 0), e(0, 1), e(0, 2));
	Vector3d X1(X(0, 0), X(0, 1), X(0, 2));
	return (e1.dot(X1));
}



MatrixXd divMatrix(std::map<int, vector<int>>& adj, MatrixXd& V, MatrixXi& F, MatrixXd& areaMatrix, MatrixXd& massMatrix, MatrixXd& gradU) {
	MatrixXd div = MatrixXd::Zero(V.rows(), 1);
	MatrixXd angles = internal_angles(V, F);
	MatrixXd X(1, 3), e1(1, 3), e2(1, 3);
	vector<int> e, ind;

	for (int i = 0; i < V.rows(); i++) {
		for (const auto& value : adj.at(i)) {
			X = gradU.row(value);
			for (int j = 0; j < 3; j = j + 1) {
				if (F(value, j) != i) {
					e.push_back(F(value, j));
					ind.push_back(j);
				}
			}
			e1 = V.row(e.at(0)) - V.row(i);
			e2 = V.row(e.at(1)) - V.row(i);
			div(i, 0) += 0.5 * cotan(angles(value, ind.at(0))) * dot(e1, X) + 0.5 * cotan(angles(value, ind.at(1))) * dot(e2, X);
			e.clear();
			ind.clear();
		}
	}


	return div;
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




int lookForNext(int sourceID, MatrixXd V, MatrixXi F, MatrixXd minCurvD, HalfedgeDS he) {
	vector<int> neighbours = vertexNeighbour(he, sourceID);

	return 0;
}

/*
int main()
{
	// Load a mesh in OFF format
	igl::opengl::glfw::Viewer viewer;
	Eigen::MatrixXd V, C;
	Eigen::MatrixXi F;

	Eigen::MatrixXd I, A;

	std::string filename = "C:\\Users\\yassi\\Desktop\\Projects\\INF574\\Geodesics_on_triangular_mesh\\data\\mctest2.off";
	igl::read_triangle_mesh(filename, V, F);


	viewer.data().set_mesh(V, F);

	HalfedgeBuilder* builder = new HalfedgeBuilder();

	HalfedgeDS he = builder->createMesh(V.rows(), F);


	VectorXi sources = VectorXi::Zero(he.sizeOfVertices());
	sources(0) = 1;
	sources(5) = 1;
	sources(10) = 1;
	sources(50) = 1;


	MatrixXd phi = dijkstra(V, F, sources);

	viewer.data().set_data(phi);
	viewer.data().show_lines = true;
	viewer.launch();

}*/


float evaluation(MatrixXd& phi, MatrixXd& realValue) {

	for (int i = 0; i < phi.rows(); i++) {
		phi(i, 0) /= phi.maxCoeff();
		realValue(i, 0) /= realValue.maxCoeff();
	}

	double error = 0;

	MatrixXd percentError = MatrixXd::Zero(phi.rows(), 1);

	for (int i = 0; i < phi.rows(); i++) {
		percentError(i, 0) = abs(phi(i, 0) - realValue(i, 0)) / realValue(i, 0);
	}
	error = percentError.sum() / phi.rows();
	error *= 100;
	return error;
}


void heatMethod(MatrixXd& V, MatrixXi& F, HalfedgeDS he) {
	const double h = igl::avg_edge_length(V, F);
	double t = pow(h, 2);
	Eigen::MatrixXd A;

	std::map<int, vector<int>> adj;



	std::cout << "calculating neighbours" << std::endl;
	buildNeighbours(adj, he, V);

	Eigen::MatrixXd U0 = Eigen::MatrixXd::Zero(V.rows(), 1);
	Eigen::MatrixXd U1 = Eigen::MatrixXd::Zero(V.rows(), 1);
	U0(0, 0) = 1;
	auto start2 = std::chrono::high_resolution_clock::now();
	auto finish2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed2 = finish2 - start2;
	std::cout << "calculation time of dijkstra: " << elapsed2.count() << " s\n";
	Eigen::SparseMatrix<double> L, M;
	Eigen::SparseMatrix<double> G;
	// Compute gradient operator:


	cout << "Calculating geodesic paths" << endl;
	igl::cotmatrix(V, F, L);
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
	igl::grad(V, F, G);

	MatrixXd div;
	A = M - t * L;
	Eigen::FullPivLU<Eigen::MatrixXd> dec(A);

	auto start = std::chrono::high_resolution_clock::now();

	U1 = dec.solve(U0);
	Eigen::MatrixXd GU = Eigen::Map<const Eigen::MatrixXd>((G * U1).eval().data(), F.rows(), 3);

	std::cout << "calculating areas" << std::endl;
	MatrixXd area;
	area = areaMatrix(V, F);

	GU = -GU;
	GU.rowwise().normalize();


	div = divMatrix(adj, V, F, area, MatrixXd(M), GU);
	Eigen::FullPivLU<Eigen::MatrixXd> dec1(L);


	phi = dec1.solve(div);
	double shift = phi.minCoeff();

	for (int i = 0; i < phi.rows(); i++) {
		phi(i, 0) -= shift;
	}
	// for measuring time performances
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "calculation time of heat method: " << elapsed.count() << " s\n";
	
}

int main()
{
	// Load a mesh in OFF format
	igl::opengl::glfw::Viewer viewer;
	Eigen::MatrixXd V, C;
	Eigen::MatrixXi F;



	std::string filename = "C:\\Users\\rached\\Documents\\Telecom\\igd\\Xinf574\\project\\data\\sphere_meshlab.off";

	igl::read_triangle_mesh(filename, V, F);
	viewer.data().set_mesh(V, F);
	
	HalfedgeBuilder* builder = new HalfedgeBuilder();
	HalfedgeDS he = builder->createMesh(V.rows(), F);



	heatMethod(V,F,he);
	
	VectorXi sources = VectorXi::Zero(he.sizeOfVertices());
	sources(0) = 1;
	phi2 = dijkstra(V, F, sources);
	
	





	viewer.callback_mouse_down = [&V, &F, &C, &he](igl::opengl::glfw::Viewer& viewer, int, int)->bool
	{
		int fid;
		Eigen::Vector3f bc;
		// Cast a ray in the view direction starting from the mouse position
		double x = viewer.current_mouse_x;


		double y = viewer.core().viewport(3) - viewer.current_mouse_y;

		if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core().view,
			viewer.core().proj, viewer.core().viewport, V, F, fid, bc))
		{

			std::cout << F(fid, 0) << " " << F(fid, 1) << " " << F(fid, 2) << std::endl;
			std::cout <<"heat method values" <<phi(F(fid, 0), 0) << " " << phi(F(fid, 1), 0) << " " << phi(F(fid, 2), 0) << std::endl;
			std::cout << "dijkstra method's values"<<phi2(F(fid, 0), 0) << " " << phi2(F(fid, 1), 0) << " " << phi2(F(fid, 2), 0) << std::endl;
			std::cout << "////////////////////////////////////////////" << std::endl;
			return true;
		}
		return false;
	};
	
	MatrixXd realValue = MatrixXd::Zero(V.rows(), 1);
	for (int i = 0; i < V.rows(); i++) { //calculating real distance for a sphere mesh
		realValue(i, 0) = 2 * acos(dot(V.row(i), V.row(0)) / 4);
	}


	std::cout << "press 1 key to visualize the heat method result" << std::endl;
	std::cout << "press 2 key to visualize the dijkstra method result" << std::endl;

	viewer.data().set_data(phi2);



	viewer.callback_key_down = &key_down;
	viewer.data().show_lines = true;
	viewer.launch();

}