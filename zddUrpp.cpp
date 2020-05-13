// general preamble
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <bitset>
#include <map>
#include <chrono>
// TdZdd-specific preamble
#include "DdSpecOp.hpp"
#include "DdStructure.hpp"
#include "util/Graph.hpp"
#include "util/IntSubset.hpp"
#include "spec/DegreeConstraint.hpp"
#include "spec/FrontierBasedSearch.hpp"

#define MAX_EDGES 200
#define MAX_REQUIRED 80

using namespace std;

// number of edges and required edges
int n=0, R=0;
// required id and weight assignments
vector<int> r, d;

// variables for MinDist optimization
vector<int> vp, vl;
vector<bool> vt;
int dd_count = 0;

// stores the data to represent the results of the ZDD
struct DistData {
	int val;
	int id;
	
	DistData(){
		val = 0;
		id = -1;
	}
	
	DistData(int v){
		val = v;
	}
	
	void print(tdzdd::Graph graph){
		cout << "Distance = " << val << endl;
		cout << "Yielding Path = ";
		
		bool first = true;
		int cur = id;
		while(cur > 0){
			if(vt[cur]){
				if(!first){
					cout << ", ";
				}else{
					first = false;
				}
				auto edgeInfo = graph.edgeInfo(n-vl[cur]);
				cout << "(" << graph.vertexName(edgeInfo.v1) << ", " << graph.vertexName(edgeInfo.v2) << ")";
			}
			
			cur = vp[cur];
		}
		
		cout << endl;
	}
};

class MinDist: public tdzdd::DdEval<MinDist,DistData> {
public:
	MinDist(){}
	void evalTerminal(DistData& data, bool one) const {
		// initialize
		data.val = one ? 0 : 100000;
	}
	void evalNode(DistData& data, int level, tdzdd::DdValues<DistData,2> const& values) const {
		// principal declaration
		vl.push_back(level);
		data.id = dd_count++;
		DistData data0 = values.get(0);
		DistData data1 = values.get(1);
		// general logic for finding the minimum
		if(data0.val <= data1.val + d[n-level]){
			data.val = data0.val;
			vt.push_back(false);
			vp.push_back(data0.id);
		}else{
			data.val = data1.val + d[n-level];
			vt.push_back(true);
			vp.push_back(data1.id);
		}
	}
};

// setting up the constraint of using each required edge at least once
class Required: public tdzdd::DdSpec<Required, int, 2> {
public:
  	Required(){}
  	int getRoot(int& state) const{
    	state = 1;
    	return n;
  	}
  	int getChild(int& state, int level, int value) const{
		if(level == n || r[n-level] != r[n-level-1]){
			if(state == 0) return 0;
			state = 0;
		}
		
		if(value == 1) state = 1;
		level--;

		if(level == 0) return -1;
    	return level;
  	}
};

int main(int argc, char *argv[]) {
	auto start_time = chrono::high_resolution_clock::now();
	
	if (argc <= 1) {
        cout << "usage: " << argv[0] << " <graph_filename>" << endl;
        return 0;
    }
	// to display computation time
	tdzdd::MessageHandler::showMessages();
	
    tdzdd::Graph graph;
	
	// read graph file
	string line, tmp;
	ifstream file;
	file.open(string(argv[1]));
	
	// trash "number of vertices : ## "
	getline(file, line);
	
	// read "number of required edges  # "
	getline(file, line);
	stringstream ss(line);
	for(int i=0; i<6; i++) getline(ss, tmp, ' ');
	
	// read required edges
	for(int i=stoi(tmp); i>0; i--){
		getline(file, line);
		stringstream ss2(line);
		
		string u, v, l;
		getline(ss2, u, ' ');
		getline(ss2, v, ' ');
		getline(ss2, l, ' ');
		
		n++;
		graph.addEdge(u, v);
		d.push_back(stoi(l));
		r.push_back(R);
		
		n++;
		graph.addEdge(v, u);
		d.push_back(stoi(l));
		r.push_back(R);
		
		R++;
	}
	
	// read "number of non required edges  # "
	getline(file, line);
	ss = stringstream(line);
	for(int i=0; i<7; i++) getline(ss, tmp, ' ');
	
	// read non required edges
	for(int i=stoi(tmp); i>0; i--){
		getline(file, line);
		stringstream ss2(line);
		
		string u, v, l;
		getline(ss2, u, ' ');
		getline(ss2, v, ' ');
		getline(ss2, l, ' ');
		
		n++;
		graph.addEdge(u, v);
		d.push_back(stoi(l));
		r.push_back(-1);
		
		n++;
		graph.addEdge(v, u);
		d.push_back(stoi(l));
		r.push_back(-1);
	}
	
	cout << "Edges: " << (n >> 1) << endl;
	cout << "Required: " << R << endl;

	file.close();
	graph.update();

    // the diagram representing all paths
    tdzdd::FrontierBasedSearch fbs(graph, 1, false, false);
	
	tdzdd::IntRange zeroOrTwo(0, 6, 2);
	// setting up degree constraint
	tdzdd::DegreeConstraint dc(graph);

	for (int v = 1; v <= graph.vertexSize(); ++v) {
		dc.setConstraint(v, &zeroOrTwo);
	}
	
	cout << "Nodes: " << graph.vertexSize() << endl;
	
	// the ZDD for all paths in general
	tdzdd::DdStructure<2> dd1(tdzdd::zddIntersection(dc, fbs));

	// the ZDD for all feasible paths
	tdzdd::DdStructure<2> dd(tdzdd::zddIntersection(dd1, Required()));

	cout << "Number of ZDD nodes = " << dd.size() << endl;
	cout << "Number of elements = " << dd.evaluate(tdzdd::ZddCardinality<>()) << endl;
	
	// minimum distance output
	MinDist dist = MinDist();
	DistData ans = dd.evaluate(dist);
	cout << "Minimum Distance = " << ans.val << endl;
	
	cout << endl << endl;
	cout << "Minimum distance tour:" << endl;
	ans.print(graph);

	auto end_time = chrono::high_resolution_clock::now();
    double time_taken = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count(); 
  
    time_taken *= 1e-9;
    cout << "Total time taken: " << time_taken << setprecision(9) << "s" << endl;
	
    return 0;
}