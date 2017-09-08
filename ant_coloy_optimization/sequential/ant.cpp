#include <bits/stdc++.h>

using namespace std;

double Pr = 0.5;
int total_ants = 11;
double cities;
double alpha = 1;
double beta = 1;
double evap_coeff = 0.5;
double Q = 50;
double min_tsp_dist = INT_MAX;
int initial_trace = 5;

vector<vector<double> >trace;
vector<vector<double> >delta_trace;
vector<vector<double> >G;
vector<vector<double> >visibility;

class ants{
	public:
	bool *visited;
	double tour_length;
	vector<double> visited_nodes;
	int current_city;

	ants(int n, int start){
		this->visited = new bool[n];
		memset(this->visited, 0, sizeof(bool)*n);
		this->tour_length = 0;
		this->visited_nodes.clear();
		
		this->current_city = start;
		this->visited_nodes.push_back(start);
		this->visited[start] = true;

	}

	void moveant_next_city(){
		double sum = 0;
		for(int s=0;s<cities;s++){
			if(!visited[s])sum = sum + pow(trace[current_city][s],alpha)*pow(visibility[current_city][s],beta);
		}
		vector<double> p;

		for(int j=0;j<cities;j++)
			if(!visited[j])p.push_back(pow(trace[current_city][j],alpha)*pow(visibility[current_city][j],beta)/sum);
			else p.push_back(-1);

		double selected_city = 0;
		vector<double> selected_city_list;
		for(int i=0;i<cities;i++)
			if(p[i]>p[selected_city]) {
				selected_city=i;
				selected_city_list.clear();
				selected_city_list.push_back(i);
			}
			else if(p[i]==p[selected_city])selected_city_list.push_back(i);

		
		visited_nodes.push_back(selected_city);
		tour_length = tour_length + G[current_city][selected_city];
		current_city = selected_city;
		visited[current_city] = true;
	}


	void update_trail(){
		double update_trace = Q/tour_length;
		for(int i=cities-1; i>0;i--)
			delta_trace[i][i-1] = delta_trace[i][i-1] + update_trace;
	}

	void update_tsp_dist(){
		if(tour_length < min_tsp_dist) min_tsp_dist = tour_length;
	}
};

vector<ants> antobj;

void setup_ants(){
	antobj.clear();
	for(int i=0;i<total_ants;i++){
		ants antobject(cities, rand()%(int)cities);
		antobj.push_back(antobject);
	}
}
	
double euclidean_dist(double x1, double y1, double x2, double y2){
	return sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}

void readGraph(char* filename){
	FILE *fp = fopen(filename,"r");
	char line[10240];
	vector<pair<double,double> >graph;

	for(int i=0;i<5;i++)fgets(line,10240,fp); //read unneccesary lines
	while(fgets(line,10240,fp)){
		double x,y,t;
		if(strcmp(line,"EOF")==0)break;
		sscanf(line,"%lf %lf %lf",&t,&x,&y);
		graph.push_back(make_pair(x,y));
	}

	cities = graph.size();
	G.resize(cities,vector<double>(cities,0));
	visibility.resize(cities, vector<double>(cities,0));
	for(int i=0;i<cities;i++)
		for(int j=i+1;j<cities;j++){
			G[i][j] = G[j][i] = euclidean_dist(graph[i].first,graph[i].second, graph[j].first, graph[j].second);
			if(G[i][j]!=0)
			visibility[i][j] = visibility[j][i] = 1/G[i][j];
		}
}

void moveant(){
	for(int j=1;j<cities;j++)
		for(int i=0;i<total_ants;i++)
			antobj[i].moveant_next_city();
}

void updatetrail(){
	for(int i=0;i<total_ants;i++)
		antobj[i].update_trail();

	for(int i=0;i<cities;i++)
		for(int j=0;j<cities;j++)
			trace[i][j] = (1-evap_coeff)*trace[i][j] + delta_trace[i][j];
}	

void update_distance(){
	for(int i=0;i<total_ants;i++)
		antobj[i].update_tsp_dist();
}		


int main(int argc, char** argv)
{
	srand(time(NULL));
	if(argc<2)return 0;
	if(argc>2)sscanf(argv[2],"%d",&total_ants);
	readGraph(argv[1]);
	trace.resize(cities, vector<double>(cities,initial_trace));
	delta_trace.resize(cities, vector<double>(cities,0));
	int q = 10;
	while(q--){
		setup_ants();
		moveant();
		for(int i=0;i<cities;i++) for(int j=0;j<cities;j++) delta_trace[i][j] = 0;
		updatetrail();
		update_distance();
		cout<<min_tsp_dist<<endl;
	}
	
	cout<<"Total Cities"<<cities<<endl;
	cout<<min_tsp_dist<<endl;
	
	return 0;
}
