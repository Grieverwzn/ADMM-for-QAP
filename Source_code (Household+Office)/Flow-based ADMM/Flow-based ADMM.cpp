// FlowbasedADMM.cpp : Defines the entry point for the console application.
//@ Jiangtao Liu, Arizona State University

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <list>
#include <omp.h>
#include "CSVParser.h"
#include <string> 
#include <algorithm>
#include <iostream>
#include <fstream>
#include <time.h>
#include <functional>
#include <stdio.h>   
#include <tchar.h>
#include <windows.h>
#include <vector>
#include <random>
#include <chrono>


using namespace std;

FILE* g_pFileDebugLog = NULL;


#define _MAX_NUMBER_OF_LINKS 5000
#define _MAX_NUMBER_OF_AGENTS 5000
#define _MAX_NUMBER_OF_ITERATIONS 1500
#define _MAX_NUMBER_OF_AGENT_TYPES 3
#define _MAX_NUMBER_OF_LINK_TYPES 4
#define _MAX_NUMBER_OF_PHYSICAL_NODES 400
#define _MAX_LABEL_COST 9999999999999


int g_number_of_agents = 0;
int g_number_of_arcs = 0;
int g_number_of_pax_groups = 0;
int g_number_of_links = 0;
int g_number_of_nodes = 0;

int g_iteration_number = 100;

using std::random_shuffle;

std::vector<int> g_pax_group_id_vector;
std::vector<float> g_group_demand_vector;
std::vector<string> g_arc_vector;
std::vector<float> g_arc_cap_vector;

float g_optimal_solution[_MAX_NUMBER_OF_ITERATIONS];
float g_link_based_agent_type_constraint_matrix[_MAX_NUMBER_OF_LINK_TYPES][_MAX_NUMBER_OF_AGENT_TYPES];
float g_link_based_agent_type_constraint_RHS[_MAX_NUMBER_OF_LINK_TYPES];

float g_agent_link_cost_matrix[_MAX_NUMBER_OF_AGENTS][_MAX_NUMBER_OF_LINKS];  // to be extended to dynamic memory allocation 
float g_agent_link_fixed_cost_matrix[_MAX_NUMBER_OF_AGENTS][_MAX_NUMBER_OF_LINKS];  // to be extended to dynamic memory allocation 
float g_agent_link_LB_cost_matrix[_MAX_NUMBER_OF_AGENTS][_MAX_NUMBER_OF_LINKS];  // to be extended to dynamic memory allocation 
float g_agent_link_UB_cost_matrix[_MAX_NUMBER_OF_AGENTS][_MAX_NUMBER_OF_LINKS];  // to be extended to dynamic memory allocation 
float rho[_MAX_NUMBER_OF_LINKS];
float subgradient_link[_MAX_NUMBER_OF_LINKS];
float lambda_link[_MAX_NUMBER_OF_LINKS];
float penalty_term[_MAX_NUMBER_OF_LINKS];//Xin
float LHS_link[_MAX_NUMBER_OF_LINKS];
float RHS_link[_MAX_NUMBER_OF_LINKS];
float LR_term[_MAX_NUMBER_OF_LINKS];
float Constraint_new_iter[_MAX_NUMBER_OF_LINKS];
float Constraint_old_iter[_MAX_NUMBER_OF_LINKS];
float GAP= 0.0;
float Best_LB = -999999999999;
float current_best_LB= -999999999999;
float Best_UB = 99999999999999;
float Best_LB_LR = -9999999999999;
float Best_UB_LR = 999999999999999;
float total_LR_B = 0;
float total_LB_1 = 0;
float total_LB = 0;
float total_UB = 0;
float max_fixed_cost = 0;
float max_rho = 0;
int** delta_pax = NULL;
int** delta_arc = NULL;
float Best_Solution = 999999999999999;
float LR_iter = 0.3;
float Shuf_iter = 0.9;
// the following data are used in agent-based ADMM

std::map<int, int> g_internal_node_seq_no_map;  // hush table, map external node number to internal node sequence no. 
std::map<int, int> g_internal_node_seq_no_to_node_id_map;  // hush table, map external node number to internal node sequence no. 
std::map<int, int> g_internal_link_no_map;  // hash table Xin
std::map<int, int> g_external_link_id_map;  // hash table

std::map<int, int> g_internal_agent_no_map;  // hash table for agents

class CLink
{
public:

	int external_link_id;
	int acting_external_link_id;
	int link_seq_no;  // internal seq no
	int from_node_seq_no;
	int to_node_seq_no;
	float fixed_cost;
	int type;
	//bool agent_type_code[_MAX_NUMBER_OF_AGENT_TYPES];
	//float agent_type_TTcost[_MAX_NUMBER_OF_AGENT_TYPES];
};


class CNode
{
public:

	int node_seq_no;  // sequence number 
	int external_node_id;      //external node number 
	double x;
	double y;
	std::vector<CLink> m_outgoing_link_vector;
};

std::vector<CNode> g_node_vector;
std::vector<CLink> g_link_vector;

class CAgent
{
public:
	int agent_id;
	int agent_type;
	int column_cost;
	float column_flow;
	int fixed_path_flag;
	string served_pax_group;
	std::vector<int> served_pax_group_vector;
	std::vector<int> set_of_allowed_links_LR; // Xin 
	std::vector<int> m_set_of_allowed_links_flag_LR; // Xin

	//string column_node_seq;
	//string column_node_time_seq;
	//std::vector<int> column_node_vector;
	//std::vector<int> column_node_time_vector;
	//std::vector<int> column_node_arc_vector;
	//std::vector<string> column_node_time_arc_vector;


	int origin_node_id;
	int destination_node_id;
	std::vector<int> path_link_seq_no_vector;
	std::vector<int> path_node_seq_no_vector;
	float path_cost;
	int Link_type_based_customized_cost_type;
	int Link_type_based_customized_cost_value;

	CAgent()
	{ 
		agent_id = 0;
		agent_type = 0;
		column_cost = 0.0;
		column_flow = 0.0;
		fixed_path_flag = 0;
		path_cost = 0;
		Link_type_based_customized_cost_type = -1;  // default value
		Link_type_based_customized_cost_value = 0;
	}
};

std::vector <CAgent> g_agent_vector;
std::vector <CAgent> g_agent_vector_best;
std::vector <int> agent_list;


class NetworkForSP  // mainly for shortest path calculation
{
public:
	int m_threadNo;  // internal thread number 

	std::list<int>  m_SENodeList;  //scan eligible list as part of label correcting algorithm 

	float* m_node_label_cost;  // label cost 
	int* m_node_predecessor;  // predecessor for nodes
	int* m_node_status_array; // update status 
	int* m_link_predecessor;  // predecessor for this node points to the previous link that updates its label cost (as part of optimality condition) (for easy referencing)
//	FILE* g_pFileDebugLog;  // file output

	float* m_link_cost_array; // link cost 
    
	std::vector<int>  m_agent_vector; // assigned agents for computing 
	std::vector<int>  m_node_vector; // assigned nodes for computing 

	NetworkForSP()
	{
	}

	void AllocateMemory(int number_of_nodes, int number_of_links)
	{
		m_node_predecessor = new int[number_of_nodes];
		m_node_status_array = new int[number_of_nodes];
		m_node_label_cost = new float[number_of_nodes];
		m_link_predecessor = new int[number_of_nodes];   // note that, the size is still the number of nodes, as each node has only one link predecessor
		m_link_cost_array = new float[number_of_links];

		for (int l = 0; l < number_of_links; l++)
		{
			m_link_cost_array[l] = g_link_vector[l].fixed_cost; //default value
		}

	}

	~NetworkForSP()
	{

		if (m_node_label_cost != NULL)
			delete m_node_label_cost;

		if (m_node_predecessor != NULL)
			delete m_node_predecessor;

		if (m_node_status_array != NULL)
			delete m_node_status_array;

		if (m_link_predecessor != NULL)
			delete m_link_predecessor;

		if (m_link_cost_array != NULL)
			delete m_link_cost_array;

	}

	// SEList: scan eligible List implementation: the reason for not using STL-like template is to avoid overhead associated pointer allocation/deallocation
	void SEList_clear()
	{
		m_SENodeList.clear();
	}

	void SEList_push_front(int node)
	{
		m_SENodeList.push_front(node);
	}

	void SEList_push_back(int node)
	{
		m_SENodeList.push_back(node);
	}

	bool SEList_empty()
	{
		return m_SENodeList.empty();
	}

	int SEList_front()
	{
		return m_SENodeList.front();
	}

	void SEList_pop_front()
	{
		m_SENodeList.pop_front();
	}

	int optimal_label_correcting(int origin_node, int destination_node, bool bUsedModifiedAgentCost_flag, int agent_a)
		// time-dependent label correcting algorithm with double queue implementation
	{
		int internal_debug_flag = 0;

		if (g_node_vector[origin_node].m_outgoing_link_vector.size() == 0)
		{
			return 0;
		}

		for (int i = 0; i < g_number_of_nodes; i++) //Initialization for all nodes
		{
			m_node_status_array[i] = 0;  // not scanned
			m_node_label_cost[i] = _MAX_LABEL_COST;
			m_node_predecessor[i] = -1;  // pointer to previous NODE INDEX from the current label at current node and time
			m_link_predecessor[i] = -1;  // pointer to previous NODE INDEX from the current label at current node and time
		}

		//Initialization for origin node at the preferred departure time, at departure time, cost = 0, otherwise, the g_A2R_simu_interval at origin node

		m_node_label_cost[origin_node] = 0;


		SEList_clear();
		SEList_push_back(origin_node);

		while (!SEList_empty())
		{
			int from_node = SEList_front();//pop a node FromID for scanning
			//cout << from_node << endl;
			SEList_pop_front();  // remove current node FromID from the SE list

			if (internal_debug_flag)
				fprintf(g_pFileDebugLog, "SP: SE node: %d\n", g_node_vector[from_node].external_node_id);

			//scan all outbound nodes of the current node
			for (int i = 0; i < g_node_vector[from_node].m_outgoing_link_vector.size(); i++)  // for each link (i,j) belong A(i)
			{
				int to_node = g_node_vector[from_node].m_outgoing_link_vector[i].to_node_seq_no;

				ASSERT(to_node <= g_number_of_nodes);
				bool  b_node_updated = false;

				int link_seq_no = g_node_vector[from_node].m_outgoing_link_vector[i].link_seq_no;

				//if (g_link_vector[link_seq_no].agent_type_code[agent_type] == false) // to be modified based on link type
				//	continue;

				//Xin 2018.12.25 
				if (g_agent_vector[agent_a].m_set_of_allowed_links_flag_LR[link_seq_no] == 0)
				{
					continue;
				}


				float link_cost = m_link_cost_array[link_seq_no];
				if (bUsedModifiedAgentCost_flag == true)
				{
					link_cost = g_agent_link_cost_matrix[agent_a][link_seq_no];  // this corresponds to point 8 in ADMM_shortest path algorithm
				}
				float new_to_node_cost = m_node_label_cost[from_node] + link_cost;


				if (internal_debug_flag)
				{
					fprintf(g_pFileDebugLog, "SP: checking from node %d, to node %d  cost = %f\n",
						g_node_vector[from_node].external_node_id,
						g_node_vector[to_node].external_node_id,
						new_to_node_cost/*, g_node_vector[from_node].m_outgoing_link_vector[i].cost*/);
				}

				if (new_to_node_cost < m_node_label_cost[to_node]) // we only compare cost at the downstream node ToID at the new arrival time t
				{

					if (internal_debug_flag)
					{
						fprintf(g_pFileDebugLog, "SP: updating node: %d current cost: %.2f, new cost %.2f\n",
							g_node_vector[to_node].external_node_id,
							m_node_label_cost[to_node], new_to_node_cost);
					}

					// update cost label and node/time predecessor

					m_node_label_cost[to_node] = new_to_node_cost;
					m_node_predecessor[to_node] = from_node;  // pointer to previous physical NODE INDEX from the current label at current node and time
					m_link_predecessor[to_node] = g_node_vector[from_node].m_outgoing_link_vector[i].link_seq_no;  // pointer to previous physical NODE INDEX from the current label at current node and time

					b_node_updated = true;

					if (internal_debug_flag)
						fprintf(g_pFileDebugLog, "SP: add node %d into SE List\n",
							g_node_vector[to_node].external_node_id);

					SEList_push_back(to_node);
					m_node_status_array[to_node] = 1;
				}

			}
		}

		if (destination_node >= 0 && m_node_label_cost[destination_node] < _MAX_LABEL_COST)
			return 1;
		else if (destination_node == -1)
			return 1;  // one to all shortest pat
		else
			return -1;


	}

	void find_path_for_agents()
	{

		    // step 1: find shortest path if needed 
			for (int i = 0; i < m_agent_vector.size(); i++)
			{

				CAgent* p_agent = &(g_agent_vector[m_agent_vector[i]]);

				if (p_agent->fixed_path_flag == 1)
					continue;

				p_agent->path_link_seq_no_vector.clear();  // reset;
				p_agent->path_node_seq_no_vector.clear();  // reset;


			   //step 3 find the destination node
				int return_value = optimal_label_correcting(g_internal_node_seq_no_map[p_agent->origin_node_id], -1, false, g_internal_agent_no_map[p_agent->agent_id]);


				if (return_value == -1)
				{
					fprintf(g_pFileDebugLog, "agent %d with can not find destination node,\n", i);
					continue;
				}

				int current_node_seq_no;
				int current_link_seq_no;

				current_node_seq_no = g_internal_node_seq_no_map[p_agent->destination_node_id];
				p_agent->path_cost = m_node_label_cost[g_internal_node_seq_no_map[p_agent->destination_node_id]];

				while (current_node_seq_no >= 0)
				{
					if (current_node_seq_no >= 0)  // this is valid node 
					{
						current_link_seq_no = m_link_predecessor[current_node_seq_no];

						if (current_link_seq_no >= 0)
						{
							p_agent->path_link_seq_no_vector.push_back(current_link_seq_no);
						}
						
						
						if (g_pFileDebugLog != NULL)
							fprintf(g_pFileDebugLog, "%d;", g_node_vector[current_node_seq_no].external_node_id);
					
						p_agent->path_node_seq_no_vector.push_back(current_node_seq_no);

					}


					current_node_seq_no = m_node_predecessor[current_node_seq_no];


				}
				if (p_agent->fixed_path_flag != 1)
				{
					std::reverse(std::begin(p_agent->path_node_seq_no_vector),
						std::end(p_agent->path_node_seq_no_vector));

					std::reverse(std::begin(p_agent->path_link_seq_no_vector),
						std::end(p_agent->path_link_seq_no_vector));
				}


				if (g_pFileDebugLog != NULL)
					fprintf(g_pFileDebugLog, "\n");


			}


	}

	void find_path_for_single_agent(int a, bool bUseModifiedCostFlag)
	{

			CAgent* p_agent = &(g_agent_vector[a]);

			if (p_agent->fixed_path_flag == 1)
				return;

			p_agent->path_link_seq_no_vector.clear();  // reset;
			p_agent->path_node_seq_no_vector.clear();  // reset;

													   //step 3 find the destination node
			int return_value = optimal_label_correcting(g_internal_node_seq_no_map[p_agent->origin_node_id], g_internal_node_seq_no_map[p_agent->destination_node_id], bUseModifiedCostFlag, a);


			if (return_value == -1)
			{
				fprintf(g_pFileDebugLog, "agent %d with can not find destination node,\n", a);
				return;
			}

			int current_node_seq_no;
			int current_link_seq_no;

			current_node_seq_no = g_internal_node_seq_no_map[p_agent->destination_node_id];
			p_agent->path_cost = m_node_label_cost[g_internal_node_seq_no_map[p_agent->destination_node_id]];

			while (current_node_seq_no >= 0)
			{
				if (current_node_seq_no >= 0)  // this is valid node 
				{
					current_link_seq_no = m_link_predecessor[current_node_seq_no];

					if (current_link_seq_no >= 0)
					{
						p_agent->path_link_seq_no_vector.push_back(current_link_seq_no);
					}

					if (g_pFileDebugLog != NULL)
						fprintf(g_pFileDebugLog, "%d;", g_node_vector[current_node_seq_no].external_node_id);

					p_agent->path_node_seq_no_vector.push_back(current_node_seq_no);

				}


				current_node_seq_no = m_node_predecessor[current_node_seq_no];


			}
			if (p_agent->fixed_path_flag != 1)
			{
				std::reverse(std::begin(p_agent->path_node_seq_no_vector),
					std::end(p_agent->path_node_seq_no_vector));

				std::reverse(std::begin(p_agent->path_link_seq_no_vector),
					std::end(p_agent->path_link_seq_no_vector));
			}


			if (g_pFileDebugLog != NULL)
				fprintf(g_pFileDebugLog, "\n");


		


	}
};

NetworkForSP* pNetworkForSP = NULL;





//split the string by ";"
vector<string> split(const string &s, const string &seperator) {
	vector<string> result;
	typedef string::size_type string_size;
	string_size i = 0;

	while (i != s.size()) {
		int flag = 0;
		while (i != s.size() && flag == 0) {
			flag = 1;
			for (string_size x = 0; x < seperator.size(); ++x)
				if (s[i] == seperator[x]) {
					++i;
					flag = 0;
					break;
				}
		}

		flag = 0;
		string_size j = i;
		while (j != s.size() && flag == 0) {
			for (string_size x = 0; x < seperator.size(); ++x)
				if (s[j] == seperator[x]) {
					flag = 1;
					break;
				}
			if (flag == 0)
				++j;
		}
		if (i != j) {
			result.push_back(s.substr(i, j - i));
			i = j;
		}
	}
	return result;
}


void g_read_network_input_data()
{
	g_number_of_agents = 0; // initialize the agent counter to 0
	g_number_of_nodes = 0; // initialize the node counter to 0
	g_number_of_links = 0;  // initialize the link counter to 0

	int internal_node_seq_no = 0;
	double x, y;
	
	// step 1: read node file 
	CCSVParser parser;
	if (parser.OpenCSVFile("input_node.csv", true))
	{
		std::map<int, int> node_id_map;

		while (parser.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
		{
			string name;
			int node_id;

			if (parser.GetValueByFieldName("node_id", node_id) == false)
				continue;
			//std::cout << node_id << endl;
			if (g_internal_node_seq_no_map.find(node_id) != g_internal_node_seq_no_map.end())
			{
				continue; //has been defined
			}
			g_internal_node_seq_no_map[node_id] = internal_node_seq_no;
			g_internal_node_seq_no_to_node_id_map[internal_node_seq_no] = node_id;
			//std::cout << g_internal_node_seq_no_map[node_id]<<endl;

			parser.GetValueByFieldName("x", x, false);
			parser.GetValueByFieldName("y", y, false);

			CNode node;  // create a node object 

			node.external_node_id = node_id;
			node.node_seq_no = internal_node_seq_no;

			node.x = x;
			node.y = y;
			internal_node_seq_no++;

			g_node_vector.push_back(node);  // push it to the global node vector
			
			g_number_of_nodes++;
			if (g_number_of_nodes % 1000 == 0)
				cout << "reading " << g_number_of_nodes << " nodes.. " << endl;
		}

		cout << "number of nodes = " << g_number_of_nodes << endl;

		fprintf(g_pFileDebugLog, "number of nodes =,%d\n", g_number_of_nodes);
		parser.CloseCSVFile();
	}
	else
	{
		cout << "input_node.csv is not opened." << endl;
		g_ProgramStop();
	}


	// step 2: read link file 
	CCSVParser parser_link;
	if (parser_link.OpenCSVFile("input_link.csv", true))
	{
		while (parser_link.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
		{
			int from_node_id = 0;
			int to_node_id = 0;
			if (parser_link.GetValueByFieldName("from_node_id", from_node_id) == false)
				continue;
			if (parser_link.GetValueByFieldName("to_node_id", to_node_id) == false)
				continue;

			// add the to node id into the outbound (adjacent) node list

			int internal_from_node_seq_no = g_internal_node_seq_no_map[from_node_id];  // map external node number to internal node seq no. 
			int internal_to_node_seq_no = g_internal_node_seq_no_map[to_node_id];

			CLink link;  // create a link object 

			parser_link.GetValueByFieldName("link_id", link.external_link_id);
			//parser_link.GetValueByFieldName("acting_link_id", link.acting_external_link_id);

			link.from_node_seq_no = internal_from_node_seq_no;
			link.to_node_seq_no = internal_to_node_seq_no;
			link.link_seq_no = g_number_of_links;
			g_internal_link_no_map[link.external_link_id] = link.link_seq_no;// Xin

			//link.to_node_seq_no = internal_to_node_seq_no;

			parser_link.GetValueByFieldName("link_type", link.type);
			parser_link.GetValueByFieldName("fixed_cost", link.fixed_cost);// distance in QAP
			max_fixed_cost = max(max_fixed_cost, link.fixed_cost);

			g_node_vector[internal_from_node_seq_no].m_outgoing_link_vector.push_back(link);  // add this link to the corresponding node as part of outgoing node/link
			g_link_vector.push_back(link);
			g_number_of_links++;

			if (g_number_of_links % 1000 == 0)
				cout << "reading " << g_number_of_links << " links.. " << endl;
		}
		parser_link.CloseCSVFile();
	}
	else
	{
		cout << "input_link.csv is not opened." << endl;
		g_ProgramStop();
	}

	cout << "number of links = " << g_number_of_links << endl;
	fprintf(g_pFileDebugLog, "number of links =,%d\n", g_number_of_links);

	// step 3: read agent/column file 
	CCSVParser parser_agent;

	if (parser_agent.OpenCSVFile("input_agent.csv", true))
	{
		while (parser_agent.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
		{
			CAgent agent;

			int r_agent_id;
			parser_agent.GetValueByFieldName("agent_id", r_agent_id);
			agent.agent_id = r_agent_id;
			g_internal_agent_no_map[agent.agent_id]= g_number_of_agents;  // hash table Xin
			
									 //Xin 2018.12.19
			int r_agent_type;
			parser_agent.GetValueByFieldName("agent_type", r_agent_type);
			agent.agent_type = r_agent_type;


			int origin_node_id = 0;
			int destination_node_id = 0;
			parser_agent.GetValueByFieldName("origin_node_id", origin_node_id, false);

			agent.origin_node_id = origin_node_id;
			parser_agent.GetValueByFieldName("destination_node_id", destination_node_id, false);
			agent.destination_node_id = destination_node_id;

			parser_agent.GetValueByFieldName("customized_cost_link_type", agent.Link_type_based_customized_cost_type, false);
			parser_agent.GetValueByFieldName("customized_cost_link_value", agent.Link_type_based_customized_cost_value, false);// Xin

			string r_column_node_seq;
			parser_agent.GetValueByFieldName("path_node_sequence", r_column_node_seq);


			// Xin input the agent-based allowed links
			int temp_link_no = 0;
			//string set_of_allowed_links_str;
			string set_of_allowed_link_types_str;
			//vector <string> set_of_allowed_links_sub_str;
			vector <string> set_of_allowed_link_types_sub_str;
			//parser_agent.GetValueByFieldName("set_of_allowed_links", set_of_allowed_links_str);
			parser_agent.GetValueByFieldName("set_of_allowed_link_types", set_of_allowed_link_types_str);
			

			for (int l = 0; l < g_number_of_links; l++)
			{
				agent.m_set_of_allowed_links_flag_LR.push_back(0);
			}

			int allowed_link_type = 0;
			if (set_of_allowed_link_types_str != "")
			{
				set_of_allowed_link_types_sub_str = split(set_of_allowed_link_types_str, ";");

				int set_of_allowed_link_type_size = max(set_of_allowed_link_types_sub_str.size(), 1);

				for (int i = 0; i < set_of_allowed_link_type_size; i++)
				{
					allowed_link_type = stoi(set_of_allowed_link_types_sub_str[i]); //stoi: string to int
					for (int l = 0; l < g_number_of_links; l++)
					{
						if (g_link_vector[l].type == allowed_link_type)
						{
							//temp_link_no = g_internal_link_no_map[l];
							temp_link_no = l;
							agent.set_of_allowed_links_LR.push_back(temp_link_no);
							agent.m_set_of_allowed_links_flag_LR[temp_link_no] = 1;
						}
					}
				}
			}
			


			g_agent_vector.push_back(agent);
			g_number_of_agents++;
		}
		parser_agent.CloseCSVFile();
	}
	cout << "number of agents = " << g_number_of_agents << endl;
	fprintf(g_pFileDebugLog, "number of agents =,%d\n", g_number_of_agents);

	// step 4: read link_type_constraint file 
	CCSVParser parser_constraint;

	for (int link_type = 0; link_type < _MAX_NUMBER_OF_LINK_TYPES; link_type++)
	{
		g_link_based_agent_type_constraint_RHS[link_type] = 999;


		for (int agent_type = 0; agent_type < _MAX_NUMBER_OF_AGENT_TYPES; agent_type++)
			g_link_based_agent_type_constraint_matrix[link_type][agent_type] = 9999;
	}

	


	if (parser_constraint.OpenCSVFile("input_constraint.csv", true))
	{
		while (parser_constraint.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
		{
			int link_type_id = 0;
			parser_constraint.GetValueByFieldName("link_type_id", link_type_id);
			if(link_type_id>=1 && link_type_id <= _MAX_NUMBER_OF_LINK_TYPES)
			{
			for (int agent_type = 0; agent_type <_MAX_NUMBER_OF_AGENT_TYPES; agent_type++)
			{
				CString agent_type_number;
				agent_type_number.Format(_T("agent_type%d"), agent_type+1);//agent_type-->agent_type+1 Xin 2018.12.17

				float coefficient = 0;
				std::string str_number = CString2StdString(agent_type_number);
				parser_constraint.GetValueByFieldName(str_number, coefficient, false);

				g_link_based_agent_type_constraint_matrix[link_type_id-1][agent_type] = coefficient;//link_type-->link_type_id-1 Xin 2018.12.17
			}
			}
			float general_capacity;
			parser_constraint.GetValueByFieldName("RHS_constant", general_capacity);
			g_link_based_agent_type_constraint_RHS[link_type_id-1] = general_capacity;//Xin
		}
		parser_constraint.CloseCSVFile();
	}
	cout << "number of constraint = " << g_number_of_links << endl;
	fprintf(g_pFileDebugLog, "number of constraint =,%d\n", g_number_of_links);
}

void g_assign_agents_to_network()
{

	pNetworkForSP = new NetworkForSP(); // create n copies of network, each for a subset of agents to use 

	pNetworkForSP->AllocateMemory(g_number_of_nodes, g_number_of_links);

	for (int a = 0; a < g_agent_vector.size(); a++)  //assign all agents to the corresponding thread
	{
		pNetworkForSP->m_agent_vector.push_back(a);
	}

	pNetworkForSP->find_path_for_agents();
}

void g_LR_ShortestPath()
{
	fprintf(g_pFileDebugLog, "start ADMM \n");
	int k, l, a;
	float step_size = 0.05;//Xin

	//point 2: initilization
	for (l = 0; l < g_link_vector.size(); l++)
	{
		subgradient_link[l] = 0;
		lambda_link[l] = 0;
	}


	for (a = 0; a < g_agent_vector.size(); a++)
	{
		for (l = 0; l < g_link_vector.size(); l++)
			g_agent_link_cost_matrix[a][l] = g_link_vector[l].fixed_cost;  //default cost
	}

	// start of ADMM
	//point 3:  outter loop for each LRiteration
	for (k = 0; k < (g_iteration_number * 8) + 1; k++)
	{

		// point 4: calculate subgradient (LHS-RHS)
		for (l = 0; l < g_link_vector.size(); l++)
		{
			LHS_link[l] = 0;
			int link_type = g_link_vector[l].type; //Xin 2018.12.19
			RHS_link[l] = g_link_based_agent_type_constraint_RHS[link_type - 1];//Xin 2018.12.19
		}

		// point 4.1 calculate (globla) LHS for each link, accumulate all links along the path of each agent 
		for (a = 0; a < g_agent_vector.size(); a++)
		{
			int agent_type = g_agent_vector[a].agent_type;
			for (int i = 0; i < g_agent_vector[a].path_link_seq_no_vector.size(); i++)
			{
				int internal_link_seq_no = g_agent_vector[a].path_link_seq_no_vector[i];
				int acting_internal_link_seq_no = g_link_vector[internal_link_seq_no].external_link_id - 1;//Xin 2018 12.25
				int link_type = g_link_vector[acting_internal_link_seq_no].type; //Xin 2018.12.19
				LHS_link[acting_internal_link_seq_no] += 1 * g_link_based_agent_type_constraint_matrix[link_type - 1][agent_type - 1];//link_type-->link_type-1 Xin 2018.12.19
			}
		}
		//}

		// point 4.3: subgradient = max(0, LHS - RHS) 
		float total_subgradient = 1;
		if (total_subgradient != 0 && k > 3)
		{
			for (l = 0; l < g_link_vector.size(); l++)
			{
				subgradient_link[l] = LHS_link[l] - RHS_link[l];
				total_subgradient += max(0, subgradient_link[l]);
				if (g_link_vector[l].type == 1)
				{
					lambda_link[l] = max(0, lambda_link[l] + step_size * subgradient_link[l]);
				}
				if (g_link_vector[l].type == 2)
				{
					lambda_link[l] = max(0, lambda_link[l] + step_size * subgradient_link[l]);
				}
				if (g_link_vector[l].type == 3)
				{
					lambda_link[l] = max(0, lambda_link[l] + step_size * subgradient_link[l]);
				}
				if (g_link_vector[l].type == 4)
				{
					lambda_link[l] = max(0, lambda_link[l] + step_size * subgradient_link[l]);
				}
			}
		}

		cout << "total_subgradient=" << total_subgradient << endl;


		//point 5: inner loop for each agent
		for (a = 0; a < g_agent_vector.size(); a++)
		{
			int agent_type = g_agent_vector[a].agent_type;


			//point 7: calculate modified cost 
			for (l = 0; l < g_link_vector.size(); l++)
			{
				//				
				float mu = LHS_link[l] - RHS_link[l]; 
				int link_type = g_link_vector[l].type;
				float alpha = g_link_based_agent_type_constraint_matrix[link_type - 1][agent_type - 1];// Xin 2018.12.24

				float customized_cost = 0;
				if (link_type == g_agent_vector[a].Link_type_based_customized_cost_type)
				{ //point 7.1 use customized agent type based cost for the specific link type
					customized_cost = g_agent_vector[a].Link_type_based_customized_cost_value;
				}

				g_agent_link_cost_matrix[a][l] = g_link_vector[l].fixed_cost * customized_cost + lambda_link[l] * alpha;//Xin			
			}
			//point 8: call shortest path with modified cost, g_agent_link_cost_matrix[a][l]  will be used in shortest path label correcting algorithm
			pNetworkForSP->find_path_for_single_agent(a, true);
			// penalty
		}  // end of point 5: inner loop for each agent
		cout << "iter_num =" << k << endl;
		fprintf(g_pFileDebugLog, "iter_num = %d\n", k);

		total_LR_B = 0;
		total_UB = 0;
		//float LR_term = 0;
		for (a = 0; a < g_agent_vector.size(); a++)
		{
			int agent_type = g_agent_vector[a].agent_type;
			for (int i = 0; i < g_agent_vector[a].path_link_seq_no_vector.size(); i++)
			{
				int internal_link_seq_no = g_agent_vector[a].path_link_seq_no_vector[i];
				int link_type = g_link_vector[internal_link_seq_no].type;
				float customized_cost = 0;

				if (link_type == g_agent_vector[a].Link_type_based_customized_cost_type)
				{ //point 7.1 use customized agent type based cost for the specific link type
					customized_cost = g_agent_vector[a].Link_type_based_customized_cost_value;
				}
				//xin
				total_LR_B += g_link_vector[internal_link_seq_no].fixed_cost * customized_cost + lambda_link[internal_link_seq_no] * g_link_based_agent_type_constraint_matrix[link_type - 1][agent_type - 1];  //default cost
				total_UB += g_link_vector[internal_link_seq_no].fixed_cost * customized_cost;  //default cost

			}
		}
		if (k > g_iteration_number)
		{
			step_size = step_size*0.999;
		}
		else
		{
			step_size = step_size;
		}
		//step_size = step_size * 0.999;

		for (l = 0; l < g_link_vector.size(); l++)
		{
			total_LR_B = total_LR_B - lambda_link[l] * RHS_link[l];
		}

		Best_LB_LR = total_LR_B;
		Best_UB_LR = total_UB;




		GAP = Best_UB_LR / Best_LB_LR - 1;
		//cout << "total_LB  =" << total_LB << endl;
		//cout << "total_UB  =" << total_UB << endl;
		cout << "Best_LB  =" << Best_LB_LR << endl;
		//cout << "Best_UB  =" << Best_UB_LR << endl;
	}// end of point 3:  outter loop for each LR iteration

}



void g_ADMM_ShortestPath()
{
	fprintf(g_pFileDebugLog, "start ADMM \n");
	int k, l, a;	
	for (l = 0; l < g_link_vector.size() + 1; l++)
	{
		rho[l] = 0.1;
	}
	float step_size = 0.1;//Xin

	//for each link
	//	lamda(l)


	//point 2: initilization
	for (l = 0; l < g_link_vector.size(); l++)
	{	subgradient_link[l] = 0; 
		//lambda_link[l] = 0.1;
		Constraint_old_iter[l] = 1;
		Constraint_new_iter[l] = 0;
	}


	for (a = 0; a < g_agent_vector.size(); a++)
	{	for (l = 0; l < g_link_vector.size(); l++)
			g_agent_link_cost_matrix[a][l] = g_link_vector[l].fixed_cost;  //default cost
			g_agent_link_LB_cost_matrix[a][l] = 0;  //default cost
			g_agent_link_UB_cost_matrix[a][l] = 0;  //default cost
			agent_list.push_back(a);
	}



	// start of ADMM
	//point 3:  outter loop for each ADMM iteration
	for (k = 0; k < g_iteration_number + 1; k++)
	{
		
		// point 4: calculate subgradient (LHS-RHS)
		for (l = 0; l < g_link_vector.size(); l++)
		{
			LHS_link[l] = 0;
			int link_type = g_link_vector[l].type; //Xin 2018.12.19
			RHS_link[l] = g_link_based_agent_type_constraint_RHS[link_type - 1];//Xin 2018.12.19
		}

			// point 4.1 calculate (globla) LHS for each link, accumulate all links along the path of each agent 
		for (a = 0; a < g_agent_vector.size(); a++)
		{
			int agent_type = g_agent_vector[a].agent_type;
			//cout <<"agent_type="<< agent_type << endl;
			for (int i = 0; i < g_agent_vector[a].path_link_seq_no_vector.size(); i++)
			{
				int internal_link_seq_no = g_agent_vector[a].path_link_seq_no_vector[i];
					//cout << internal_link_seq_no << endl;
				int acting_internal_link_seq_no = g_link_vector[internal_link_seq_no].external_link_id-1;//Xin 2018 12.25
				int link_type = g_link_vector[acting_internal_link_seq_no].type; //Xin 2018.12.19
					//int link_type = g_link_vector[l].type;   //Xin 2018.12.19

					// point 4.2 // alpha coefficeint (link_type,agent_type)
				LHS_link[acting_internal_link_seq_no] += 1 * g_link_based_agent_type_constraint_matrix[link_type-1][agent_type-1];//link_type-->link_type-1 Xin 2018.12.19
			}
		}
		//}

			// point 4.3: subgradient = max(0, LHS - RHS) 
		float total_subgradient = 0;

			for (l = 0; l < g_link_vector.size(); l++)
			{
				subgradient_link[l] = LHS_link[l] - RHS_link[l];
				total_subgradient += max(0, subgradient_link[l]);
				if ((g_link_vector[l].type == 4)||(g_link_vector[l].type == 1))
				{
					lambda_link[l] = max(0, lambda_link[l] + rho[l] * subgradient_link[l]);
					//lambda_link[l] = max(0, lambda_link[l] + step_size * subgradient_link[l]);
				}
				if ((g_link_vector[l].type == 2) || (g_link_vector[l].type == 3))
				{
					lambda_link[l] = max(0, lambda_link[l] + rho[l] * subgradient_link[l]);
				}
			}
		//cout <<"total_subgradient="<< total_subgradient << endl;
		if (k >= (g_iteration_number + 1)*Shuf_iter)
		{
			random_shuffle(agent_list.begin(), agent_list.end());
		}
		//point 5: inner loop for each agent
		for (int a1 = 0; a1 < g_agent_vector.size(); a1++)
		{
			a = agent_list[a1];
			//cout << a << endl;
			int agent_type = g_agent_vector[a].agent_type;
			

			//point 6: calculate inner-loop (local) LHS
			for (int i = 0; i < g_agent_vector[a].path_link_seq_no_vector.size(); i++)
			{
				int internal_link_seq_no = g_agent_vector[a].path_link_seq_no_vector[i];  // from previous iteration
				int acting_internal_link_seq_no = g_link_vector[internal_link_seq_no].external_link_id - 1;//Xin 2018 12.25
				int link_type = g_link_vector[acting_internal_link_seq_no].type; //Xin 2018.12.19
				//int link_type = g_link_vector[l].type;
				// point 6.1 // alpha coefficeint (link_type,agent_type)
				LHS_link[acting_internal_link_seq_no] -= 1* g_link_based_agent_type_constraint_matrix[link_type-1][agent_type-1]; // used in ADMM innner loop to represent all the other vehicles // Xin link_type-1 2018.12.19
			}

			//point 7: calculate modified cost 
			for (l = 0; l < g_link_vector.size(); l++)
			{
				//				
				float mu = LHS_link[l] - RHS_link[l];  //please recall that, LHS_link[l] is the LHS without my own coefficient 
				int link_type = g_link_vector[l].type;
				float alpha = g_link_based_agent_type_constraint_matrix[link_type-1][agent_type-1];// Xin 2018.12.24
				
				float customized_cost = 0;
				if (link_type == g_agent_vector[a].Link_type_based_customized_cost_type)
				{ //point 7.1 use customized agent type based cost for the specific link type
					customized_cost= g_agent_vector[a].Link_type_based_customized_cost_value;
				}
				//point 7.2: all overall cost = agent-based link cost + Lambda * alpha coefficient + rho *(2 * mu*alpha + alpha*alpha)
				penalty_term[l] = 0;
				if (k >= (g_iteration_number + 1)*0)
				{
					if ((link_type == 4)||(link_type == 1))
					{
						//penalty_term[l] = 0.5*rho[l] * (2 * mu*alpha + alpha * alpha);//xin
						penalty_term[l] = 0.5*rho[l] * (2 * mu*alpha + alpha * alpha);//xin
					}
					if ((link_type == 2) || (link_type == 3))
					{
						penalty_term[l] = 0.5*rho[l] * (2 * mu*alpha + alpha * alpha);//xin
					}
				}
				g_agent_link_cost_matrix[a][l] = g_link_vector[l].fixed_cost * customized_cost + lambda_link[l] * alpha +penalty_term[l];//Xin			
			}
			//point 8: call shortest path with modified cost, g_agent_link_cost_matrix[a][l]  will be used in shortest path label correcting algorithm
			pNetworkForSP->find_path_for_single_agent(a, true);

			//point 9: update local LHS again based on the shortest path calculation result
			for (int i = 0; i < g_agent_vector[a].path_link_seq_no_vector.size(); i++)
			{
				int internal_link_seq_no = g_agent_vector[a].path_link_seq_no_vector[i];  //from shortest path
				//int link_type = g_link_vector[l].type;
				int acting_internal_link_seq_no = g_link_vector[internal_link_seq_no].external_link_id - 1;//Xin 2018 12.25
				int link_type = g_link_vector[acting_internal_link_seq_no].type; //Xin 2018.12.19
				//int link_type = g_link_vector[internal_link_seq_no].type;
				float alpha = g_link_based_agent_type_constraint_matrix[link_type-1][agent_type-1];// Xin modifies
				LHS_link[acting_internal_link_seq_no] += 1* alpha;  // used in ADMM innner loop
			}
			// penalty
		}  // end of point 5: inner loop for each agent
		cout << "iter_num =" << k << endl;
		fprintf(g_pFileDebugLog, "iter_num = %d\n", k);

		for (l = 0; l < g_link_vector.size(); l++)
		{
			Constraint_new_iter[l] = LHS_link[l] - RHS_link[l];
			float c1 = abs(Constraint_new_iter[l]* Constraint_new_iter[l]);
			float c2 = abs(Constraint_old_iter[l]* Constraint_old_iter[l]);
			float c3 = c1 / c2;
			
			if (c3 >= 0.25)//(Constraint_new_iter[l]>0 && Constraint_new_iter[l] / Constraint_old_iter[l] >= 1)
			{
				//if (k >= (g_iteration_number + 1)*LR_iter)
				//{
				//	float step_rand = rand() % 100 / (double)101;
				//	if (total_subgradient>0)
				//	{
				//		rho[l] = rho[l] + step_rand;
				//	}
				//	else{
				//		if (subgradient_link[l] <= 0.00000001)
				//		{
				//			rho[l] = rho[l]* step_rand;
				//		}
				//		if (subgradient_link[l] > 0.00000001)
				//		{
				//			rho[l] = rho[l] + 1;// step1;//0.08;
				//		}
				//	}
				//}

				if (k >= (g_iteration_number + 1)*LR_iter)
				{
					rho[l] = rho[l] * 1;
					float step_rand = rand() % 100 / (double)101;
					if (total_subgradient>0)
					{
						rho[l] = rho[l] + 10*step_rand*1;
					}
					else
					{
						rho[l] = rho[l] * step_rand;
					}
				}


			}



			Constraint_old_iter[l] = Constraint_new_iter[l];
			max_rho = max(max_rho, rho[l]);
			LR_term[l] = 0; //initialization LR term
		}

		
		
		total_LR_B = 0;
		total_UB = 0;
		
		//float LR_term = 0;
		for (a = 0; a < g_agent_vector.size(); a++)
		{

			int agent_type = g_agent_vector[a].agent_type;
			for (int i = 0; i < g_agent_vector[a].path_link_seq_no_vector.size(); i++)
			{
				int internal_link_seq_no = g_agent_vector[a].path_link_seq_no_vector[i];
				int link_type = g_link_vector[internal_link_seq_no].type;
				float customized_cost = 0;

				if (link_type == g_agent_vector[a].Link_type_based_customized_cost_type)
				{ //point 7.1 use customized agent type based cost for the specific link type
					customized_cost = g_agent_vector[a].Link_type_based_customized_cost_value;
					//cout << customized_cost << endl;
				}
				//xin
				//LR_term[internal_link_seq_no] = 0;
				LR_term[internal_link_seq_no] += lambda_link[internal_link_seq_no] * g_link_based_agent_type_constraint_matrix[link_type - 1][agent_type - 1];
				//cout << "link"<< internal_link_seq_no<<"=" << lambda_link[internal_link_seq_no] * g_link_based_agent_type_constraint_matrix[link_type - 1][agent_type - 1] << endl;
				total_LR_B += g_link_vector[internal_link_seq_no].fixed_cost * customized_cost + lambda_link[internal_link_seq_no] * g_link_based_agent_type_constraint_matrix[link_type - 1][agent_type - 1];  //default cost
				//total_LR_B += lambda_link[internal_link_seq_no] * g_link_based_agent_type_constraint_matrix[link_type - 1][agent_type - 1];  //default cost
				total_UB += g_link_vector[internal_link_seq_no].fixed_cost * customized_cost; //default cost
				
			}
			//cout << total_UB << endl;
		}

		//step_size = step_size/(k+1);
		step_size = step_size*0.999;
		total_subgradient = 0;
		for (l = 0; l < g_link_vector.size(); l++)
		{
			subgradient_link[l] = LHS_link[l] - RHS_link[l];
			total_subgradient += max(0, subgradient_link[l]);
			//cout << "subgradient=" << l << ":" << subgradient_link[l] << endl;
			//cout << "total subgradient=" << total_subgradient << endl;
		}
		cout << "total_subgradient=" << total_subgradient << endl;

		float total_LB_1 = 0;
		total_LB_1 = total_LR_B;

		for (l = 0; l < g_link_vector.size(); l++)
		{	
			total_LB_1 = total_LB_1 - lambda_link[l]* RHS_link[l];
			//cout << "LR_multiplier"<<l<<"=" << lambda_link[l] << endl;
			//cout << "term" << l << "  1=" <<LR_term[l] << endl;
			//cout << "term" << l << "  2=" <<lambda_link[l] * RHS_link[l] << endl;
			//cout << "term" << l << "  3=" <<LR_term[l]- lambda_link[l] * RHS_link[l] << endl;
		}

		total_UB = total_UB;// +total_subgradient * max_fixed_cost;
		Best_UB = total_UB;
		if (k <= (g_iteration_number + 1)*LR_iter)
		{
			if (current_best_LB<=total_LB_1)
			{
				current_best_LB = total_LB_1;
				Best_LB = current_best_LB;
			}
		}

		if (total_subgradient<=0.0000000000000000001)
		{ 
			if (Best_UB <= Best_Solution)
			{
				GAP = (Best_Solution - Best_LB)/ Best_Solution;
				Best_Solution = Best_UB;
				g_agent_vector_best = g_agent_vector;
				cout << "Best_Solution  =" << Best_Solution << endl;
			}
		}

		cout << "Best_LB  =" << Best_LB << endl;
		cout << "Best_UB  =" << Best_UB << endl;
		cout << "Best_Solution  =" << Best_Solution << endl;
		cout << "GAP =" << GAP << endl;

		fprintf(g_pFileDebugLog, "Best_LB = %0.2f\n", Best_LB);
		fprintf(g_pFileDebugLog, "Best_UB = %0.2f\n", Best_UB);
		fprintf(g_pFileDebugLog, "Best_Solution = %0.2f\n", Best_Solution);
		fprintf(g_pFileDebugLog, "total_subgradient = %0.2f\n", total_subgradient);
		fprintf(g_pFileDebugLog, "GAP = %0.4f\n", GAP);



	}// end of point 3:  outter loop for each ADMM iteration

}



void g_output_agent_solutions()
{
g_agent_vector = g_agent_vector_best; // xin 0710/2019

FILE* g_pFileAgent = NULL;
g_pFileAgent = fopen("output_agent.csv", "w");
if (g_pFileAgent == NULL)
{
	cout << "File output_agent.csv cannot be opened." << endl;
	g_ProgramStop();
}
else
{

	fprintf(g_pFileAgent, "agent_id,origin_node_id,destination_node_id,cost,path_node_sequence,path_link_sequence\n");

	for (int a = 0; a < g_agent_vector.size(); a++)
	{
		CAgent* p_agent = &(g_agent_vector[a]);

		fprintf(g_pFileAgent, "%d,%d,%d,%f,",
			p_agent->agent_id,
			p_agent->origin_node_id,
			p_agent->destination_node_id,
			p_agent->path_cost
			);

		// path node id sequence
		for (int i = 0; i < p_agent->path_node_seq_no_vector.size(); i++)
		{
			int internal_node_id = p_agent->path_node_seq_no_vector[i];
			int external_node_id = g_node_vector[internal_node_id].external_node_id;
			fprintf(g_pFileAgent, "%d;", external_node_id);
		}

		fprintf(g_pFileAgent, ",");
		// path link id sequence 
		for (int i = 0; i < p_agent->path_link_seq_no_vector.size(); i++)
		{
			int internal_link_id = p_agent->path_link_seq_no_vector[i];
			int external_link_id = g_link_vector[internal_link_id].external_link_id;

			if (external_link_id >= 1)
			{
				fprintf(g_pFileAgent, "%d;", external_link_id);
			}
		}
		fprintf(g_pFileAgent, ",");
		fprintf(g_pFileAgent, "\n");

	}
	//
	fclose(g_pFileAgent);

}
}

// Main loop
int main()
{
	g_pFileDebugLog = fopen("Debug.txt", "w");
	if (g_pFileDebugLog == NULL)
	{
		cout << "File Debug.txt cannot be opened." << endl;
		g_ProgramStop();
	}
	else
	{
		cout << "File Debug.txt has been opened" << endl;
	}
	
	cout << "reading input_data" << endl;
	// Step 1 
	g_read_network_input_data();
	// Step 2
	g_assign_agents_to_network();
	//g_LR_ShortestPath();

	g_ADMM_ShortestPath(); //Xin 2018.12.17
	g_output_agent_solutions();


	cout << "well done!" << endl;

}

