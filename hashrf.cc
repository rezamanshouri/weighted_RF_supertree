
#include <sys/time.h>
#include <sys/resource.h>
#include <cassert>
#include <valarray>
#include <fstream>
#include <iostream>

#include "label-map.hh"
#include "hashfunc.hh"
#include "hash.hh"
#include "./tclap/CmdLine.h"


//////////////////////REZA////////////////////////////////////////>>>>>>>>>>>>>>>>>>>>
//////////////////////////////////////////////////////////////////

#include <regex>
#include <iterator>
#include <set>
#include <string>

void split(const string &s, char delim, vector<string> &elems);
vector<string> split(const string &s, char delim);
set<string> find_non_common_taxa_set(const string &supertree, const string &source_tree);
//void restrict_supertree_to_source_tree(NEWICKNODE *node, const set<string> non_shared_taxa);
string change_branch_lengths(const string &tree, int percentage_to_be_reweighted, int new_weight);
void write_line_to_File(string s,const char* file_name);
double calculate_wrf_btwn_ST_n_source_tree(int argc, char** argv);
void restrict_st_to_source_tree();
double calculate_total_wrf(const char* input_file);
double calculate_total_rf(const char* input_file);


//************ADDED FOR SPR neihborhood************>>>>>>>>>>>>>>>>>>>>>>>>>>>
#include <cstdio>
#include <cstdlib>
#include <ctime>
//#include <string>
#include <cstring>
//#include <iostream>
//#include <fstream>
#include <sstream>
#include <climits>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <list>
#include <time.h>
#include "rspr.h"

#include "Forest.h"
#include "ClusterForest.h"
#include "LCA.h"
#include "ClusterInstance.h"
#include "UndoMachine.h"
#include "lgt.h"
#include "sparse_counts.h"
#include "node_glom.h"


#include <regex>
void writeToFile(string s,const char* file_name);
void adjustTree(Node* tree);
int produce_all_spr_neighbors(Node* tree, int number_of_taxa);
int produce_n_percent_of_spr_neighbors(Node* myTree, int number_of_taxa, int percentage);
int number_of_taxa(string const tree);
void total_number_of_nodes(Node* node, int& total_nodes);

////////////////////////REZA////////////////////////////////////////////<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
////////////////////////////////////////////////////////////////////        






// For newick parser
extern "C" {
#include <newick.h>
}

using namespace std;

#define LEFT               								0
#define RIGHT              								1
#define ROOT               								2
#define LABEL_WIDTH        								3

// Set a random number for m1 (= Initial size of hash table)
// m1 is the closest value to (t*n)+(t*n*HASHTABLE_FACTOR)
#define HASHTABLE_FACTOR    							0.2

// the c value of m2 > c*t*n in the paper
static unsigned int C						    			= 1000; 
static unsigned NUM_TREES                 = 0; // number of trees
static unsigned NUM_TAXA                	= 0; // number of taxa
static bool WEIGHTED                      = false; // unweighted
static unsigned PRINT_OPTIONS							= 3; // matrix
static int32 NEWSEED											= 0; // for storing user specified seed value 

#define __HALF_MAX_SIGNED(type) ((type)1 << (sizeof(type)*8-2))
#define __MAX_SIGNED(type) (__HALF_MAX_SIGNED(type) - 1 + __HALF_MAX_SIGNED(type))
#define __MIN_SIGNED(type) (-1 - __MAX_SIGNED(type))
#define __MIN(type) ((type)-1 < 1?__MIN_SIGNED(type):(type)0)
#define __MAX(type) ((type)~__MIN(type))

#define assign(dest,src) ({ \
  decltype(src) __x=(src); \
  decltype(dest) __y=__x; \
  (__x==__y && ((__x<1) == (__y<1)) ? (void)((dest)=__y),0:1); \
})

#define add_of(c,a,b) ({ \
  decltype(a) __a=a; \
  decltype(b) __b=b; \
  (__b)<1 ? \
    ((__MIN(decltype(c))-(__b)<=(__a)) ? assign(c,__a+__b):1) : \
    ((__MAX(decltype(c))-(__b)>=(__a)) ? assign(c,__a+__b):1); \
})

void
GetTaxaLabels(
  NEWICKNODE *node,
  LabelMap &lm)
{
  if (node->Nchildren == 0) {
    string temp(node->label);
    lm.push(temp);
  }
  else
    for (int i=0;i<node->Nchildren;++i)
      GetTaxaLabels(node->child[i], lm);
}

void
dfs_compute_hash(
  NEWICKNODE* startNode,
  LabelMap &lm,
  HashRFMap &vvec_hashrf,
  unsigned treeIdx,
  unsigned &numBitstr,
  unsigned long long m1,
  unsigned long long m2)
{
  // If the node is leaf node, just set the place of the taxon name in the bit string to '1'
  // and push the bit string into stack
  if (startNode->Nchildren == 0) { // leaf node
    string temp(startNode->label);
    unsigned int idx = lm[temp];

    // Implicit BPs /////////////////////////
    // Set the hash values for each leaf node.
    startNode->hv1 = vvec_hashrf._HF.getA1(idx);
    startNode->hv2 = vvec_hashrf._HF.getA2(idx);
  }
  else {
    for (int i=0; i<startNode->Nchildren; ++i) {
      dfs_compute_hash(startNode->child[i], lm, vvec_hashrf, treeIdx, numBitstr, m1, m2);
    }

    // For weighted RF
    float dist = 0.0;
    if (WEIGHTED)
      dist = startNode->weight;
    else
      dist = 1;

    ++numBitstr;

    // Implicit BPs ////////////
    // After an internal node is found, compute the hv1 and hv2
    unsigned long long temphv1=0;
    unsigned long long temphv2=0;
   
    for (int i=0; i<startNode->Nchildren; ++i) {    	
    	unsigned long long t1 = temphv1;
    	unsigned long long t2 = temphv2;
    	unsigned long long h1 = startNode->child[i]->hv1;
    	unsigned long long h2 = startNode->child[i]->hv2;

    	if ( add_of(temphv1, t1, h1) ) {
    		cout << "ERROR: ullong add overflow!!!\n"; 
    		cout << "t1=" << t1 << " h1=" << h1 << " t1+h1=" << t1+h1 << endl;
    		exit(0);
    	}
    	if ( add_of(temphv2, t2, h2) ) {
    		cout << "ERROR: ullong add overflow!!!\n"; 
    		cout << "t2=" << t2 << " h2=" << h2 << " t2+h2=" << t2+h2 << endl;
    		exit(0);
    	}
	  }
		
		// Check overflow
		unsigned long long temp1 = temphv1 % m1;
		unsigned long long temp2 = temphv2 % m2;
    	startNode->hv1 = temp1;
    	startNode->hv2 = temp2;

    // Store bitstrings in hash table
    if (numBitstr < NUM_TAXA-2) {
      vvec_hashrf.hashing_bs_without_type2_nbits(treeIdx, NUM_TAXA, startNode->hv1, startNode->hv2, dist, WEIGHTED);   // without TYPE-III using n-bits (Hash-RF)
    }
  }
}



static double
print_rf_short_matrix(
  vector< vector<unsigned short> > &SHORTSIM,
  unsigned options,
  string outfile)
{
  double rf_dist= 0.0;

  ofstream fout;
  if (outfile != "")
  {
    fout.open(outfile.c_str());
  }

  switch (options) {
  case 0:
    return -1;
  case 1:
    cout << "\nRobinson-Foulds distance (list format):\n";

    if (outfile == "") {
      for (unsigned int i = 0; i < NUM_TREES; ++i) {
        for (unsigned int j = 0; j < NUM_TREES; ++j) {
          if (i == j)
            cout << "<" << i << "," << j << "> " << 0 << endl;
          else {
              cout << "<" << i << "," << j << "> " << (NUM_TAXA-3)-(float((SHORTSIM[i][j] + SHORTSIM[j][i])/2)) << endl;
          }
        }
      }
    }
    else {
      for (unsigned int i = 0; i < NUM_TREES; ++i) {
        for (unsigned int j = 0; j < NUM_TREES; ++j) {
          if (i == j)
            fout << "<" << i << "," << j << "> " << 0 << endl;
          else {
              fout << "<" << i << "," << j << "> " << (NUM_TAXA-3)-(float((SHORTSIM[i][j] + SHORTSIM[j][i])/2)) << endl;
          }
        }
      }
    }

    break;
  case 2:
    cout << "\nRobinson-Foulds distance (rate):\n";
    if (WEIGHTED) {cout << "Fatal error: RF rate is only for unweighted RF distance.\n"; exit(0);}

    if (outfile == "") {
      for (unsigned int i = 0; i < NUM_TREES; ++i) {
        for (unsigned int j = 0; j < NUM_TREES; ++j) {
          cout << "<" << i << "," << j << "> ";
          if (i==j)
            cout << 0 << endl;
          else
            cout << (float) ((NUM_TAXA-3)-((SHORTSIM[i][j] + SHORTSIM[j][i])/2)) / (NUM_TAXA-3) * 100 << endl;
        }
      }
    }
    else {
      for (unsigned int i = 0; i < NUM_TREES; ++i) {
        for (unsigned int j = 0; j < NUM_TREES; ++j) {
          fout << "<" << i << "," << j << "> ";
          if (i==j)
            fout << 0 << endl;
          else
            fout << (float) ((NUM_TAXA-3)-((SHORTSIM[i][j] + SHORTSIM[j][i])/2)) / (NUM_TAXA-3) * 100 << endl;
        }
      }
    }
    break;
  case 3:
    //cout << "\nRobinson-Foulds distance (matrix format):\n";
    if (outfile == "") {
      rf_dist = (NUM_TAXA-3)-(float((SHORTSIM[0][1] + SHORTSIM[1][0])/2));
      //cout << "non_w: " << rf_dist << endl;
      return rf_dist;
      /*
      for (unsigned int i = 0; i < NUM_TREES; ++i)  {
        for (unsigned int j = 0; j < NUM_TREES; ++j)  {
//	        for (unsigned int j = i; j < NUM_TREES; ++j)  {
          if (i == j)
            cout << "0" << ' ';
          else
            cout << (NUM_TAXA-3)-(float((SHORTSIM[i][j] + SHORTSIM[j][i])/2)) << ' ';
        }
        cout << endl;
      }
      cout << endl;
      */
    }
    else {
      for (unsigned int i = 0; i < NUM_TREES; ++i)  {
        for (unsigned int j = 0; j < NUM_TREES; ++j)  {
          if (i == j)
            fout << "0" << ' ';
          else
            fout << (NUM_TAXA-3)-(float((SHORTSIM[i][j] + SHORTSIM[j][i])/2)) << ' ';
        }
        fout << endl;
      }
      fout << endl;
    }
    break;
  case 4:
  	if (outfile == "") { 
	    for (size_t i = 0; i < NUM_TREES; ++i) {		
		    for	(size_t j = 0; j < i; ++j)
		      std::cout << (NUM_TAXA-3)-(float((SHORTSIM[i][j] + SHORTSIM[j][i])/2)) << " ";
		    std::cout << "\n";
		 	}		
		}
		else {
			for (size_t i = 0; i < NUM_TREES; ++i) {		
		    for	(size_t j = 0; j < i; ++j)
		      fout << (NUM_TAXA-3)-(float((SHORTSIM[i][j] + SHORTSIM[j][i])/2)) << " ";
		    	fout << "\n";
		 	}	
		}
    break;
  }

  if (outfile != "")
    fout.close();
}
  
//I change thid function to return distance instead of printting it
//note I ALWAYS have exactly two trees: ST and one of source trees  
static double
print_rf_float_matrix(
	vector< vector<float> > &SIM, 
	unsigned options,
	string outfile)
{
  double wrf_dist= 0.0;

	ofstream fout;
	if (outfile != "") {
		fout.open(outfile.c_str());
	}
	
	switch (options) {
		case 0: 
			//return;
      return -1;
      cout << "something went wrong" << endl;   
		case 1: 
			cout << "\nRobinson-Foulds distance (list format):\n";
			
				if (outfile == "") {
					for (unsigned i = 0; i < NUM_TREES; ++i) {          
						for (unsigned j = 0; j < NUM_TREES; ++j) {        	
							if (i == j) 
								cout << "<" << i << "," << j << "> " << 0 << endl;
							else {
								if (!WEIGHTED)
									cout << "<" << i << "," << j << "> " << (NUM_TAXA-3)-(float((SIM[i][j] + SIM[j][i])/2)) << endl;
								else 
									cout << "<" << i << "," << j << "> " << float((SIM[i][j] + SIM[j][i])/4) << endl;
							}
						}
					}
				}
				else {
					for (unsigned i = 0; i < NUM_TREES; ++i) {          
						for (unsigned j = 0; j < NUM_TREES; ++j) {        	
							if (i == j) 
								fout << "<" << i << "," << j << "> " << 0 << endl;
							else {
								if (!WEIGHTED) 
									fout << "<" << i << "," << j << "> " << (NUM_TAXA-3)-(float((SIM[i][j] + SIM[j][i])/2)) << endl;
								else 
									fout << "<" << i << "," << j << "> " << float((SIM[i][j] + SIM[j][i])/4) << endl;
							}
						}
					}
				}
			
			break;      
		case 2:
			cout << "\nRobinson-Foulds distance (rate):\n"; 
			if (WEIGHTED) {
					cout << "Fatal error: RF rate is only for unweighted RF distance.\n";
					exit(0);
			}      
			
			if (outfile == "") {
				for (unsigned i = 0; i < NUM_TREES; ++i) {
					for (unsigned j = 0; j < NUM_TREES; ++j) {              
						cout << "<" << i << "," << j << "> ";
						if (i==j) 
							cout << 0 << endl;
						else 
							cout << (float) ((NUM_TAXA-3)-((SIM[i][j] + SIM[j][i])/2)) / (NUM_TAXA-3) * 100 << endl;
					}
				}
			}
			else {
				for (unsigned i = 0; i < NUM_TREES; ++i) {
					for (unsigned j = 0; j < NUM_TREES; ++j) {              
						fout << "<" << i << "," << j << "> ";
						if (i==j) 
							fout << 0 << endl;
						else      
							fout << (float) ((NUM_TAXA-3)-((SIM[i][j] + SIM[j][i])/2)) / (NUM_TAXA-3) * 100 << endl;
					}
				}
			}
			break;      
		case 3:
			//cout << "\nRobinson-Foulds distance (matrix format):\n";
			if (outfile == "") {
         wrf_dist =  float((SIM[0][1] + SIM[1][0])/4);
         //cout << "w_dits " << wrf_dist << endl;
         return wrf_dist;
        /*
				for (unsigned i = 0; i < NUM_TREES; ++i)  {
					for (unsigned j = 0; j < NUM_TREES; ++j)  {
//	        for (unsigned j = i; j < NUM_TREES; ++j)  {
						if (i == j)
							cout << " 0 " << ' ';
						else
							if (WEIGHTED)
								cout << float((SIM[i][j] + SIM[j][i])/4) << ' ';
							else 
								cout << (NUM_TAXA-3)-(float((SIM[i][j] + SIM[j][i])/2)) << ' ';
					}
					cout << endl;
				}
				cout << endl;
        */
			}
			else {
				for (unsigned i = 0; i < NUM_TREES; ++i)  {
					for (unsigned j = 0; j < NUM_TREES; ++j)  {
						if (i == j)
							fout << " 0 " << ' ';
						else
							if (WEIGHTED)
								fout << float((SIM[i][j] + SIM[j][i])/4) << ' ';
							else 
								fout << (NUM_TAXA-3)-(float((SIM[i][j] + SIM[j][i])/2)) << ' ';
					}
					fout << endl;
				}
				fout << endl;
			}
			break;    
	}
	
	if (outfile != "") 
		fout.close();	
}









//////////////////////////////////////REZA///////////////////////////////main()
//////////////////////////////////////////////////////>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
int main(int argc, char** argv) {
  srand((unsigned)time(NULL));

  //Remember each ratchtet iteration has two separate hill-climbing seraches
  //1- on weighted characters: for this step, I use "hashrf z_weighted_input_trees_reweighted 0 -w" command, where z_weighted_input_trees file is generated below and modified in each iteration 
  //2- on unweighted characters: for this step, I use "hashrf argv[1] 0" command
  
  //Note for varaiable names, I used "..._wrf_dist_..." which is sometimes actually wrf_dist and sometimes rf_dist
  //based on whether we are in what hill_climbing phase of each ratchet iter. This is fine since
  //the way both are done is the same. However, the_best_supertree_seen is updated only when we are in
  //unweighted hill-climbing phase of the ratchet iter and found a better rf-dist.


  //file "z_weighted_input_trees_reweighted" is created from "z_weighted_input_trees_original" by calling 
  // change_branch_lengths() on each line of it and copying the output to "z_weighted_input_trees_reweighted"

  //add weight 1 to all internal edges: using perl replace all ")" with ""):1" except last one i.e. "):1;"
  string argv1(argv[1]);
  string command1 = "perl -pe 's/\\)/):1/g' " + argv1 + " > z_weighted_input_trees_original";
  system(command1.c_str());
  string command2 = "perl -pe 's/:1;/;/g' z_weighted_input_trees_original > z_temp; cat z_temp> z_weighted_input_trees_original; rm z_temp";
  system(command2.c_str());
  
  //Note we don't need to "weight" ST since when calling restrict_st_to_source_tree(), it gets weighted.




  cout << "Please enter the pre-specified number of 'ratchet' iterations: " << endl;
  //pre-specified number of iterations if it didn't stop after this number of iterations
  int number_of_ratchet_iterations;
  cin >> number_of_ratchet_iterations;

  cout << "Please enter the percentage of clades to be re-weighted in Ratchet search (between 0 and 100): " << endl;
  int percentage_of_clades_to_be_reweighted;
  cin >> percentage_of_clades_to_be_reweighted;

  cout << "Please enter the new weight to which clades be re-weighted in Ratchet search (between 0 and 10): " << endl;
  int ratchet_weight;
  cin >> ratchet_weight;  

  cout << "Please enter the percentage of the SPR-neighborhood to be considered (randomly) in each local search (between 0 and 100): " << endl;
  int neighborhood_percentage;
  cin >> neighborhood_percentage; 


  string init_supertree;
  ifstream init_sup("initial_supertree");
  if (init_sup.good())
  {
      string sLine;
      getline(init_sup, sLine);
      init_supertree= sLine;
  }else{
      cout << "Initial supertree file not found!!" << endl;
      return -1;
  }

  cout << "********************************************************************" << endl;



  //running time
  clock_t start_time,finish_time;
  start_time= clock();


  //just to keep track of best ST seen throughout the algorithm
  string the_best_supertree_seen = init_supertree;
  int the_best_wrf_distance_seen = INT_MAX;


  ////////////////////////////////////////////////
  //////////////////////ratchet///////////////////
  ////////////////////////////////////////////////

  //ratchet search loop
  for(int ratchet_counter = 1; ratchet_counter < number_of_ratchet_iterations+1; ratchet_counter++) {

      bool ratchet= true; //when true, use re-weighted input trees; when false, use unweighted trees.

      const char* source_trees_for_this_search;
      if(ratchet) {
        ofstream reweighted_input_trees;
        reweighted_input_trees.open("z_weighted_input_trees_reweighted");
        ifstream original_weighted_trees;
        original_weighted_trees.open("z_weighted_input_trees_original");
        string current_tree;
        while(getline(original_weighted_trees, current_tree)) {
          current_tree = change_branch_lengths(current_tree, percentage_of_clades_to_be_reweighted, ratchet_weight) + "\n";
          reweighted_input_trees << current_tree;          
        }
        original_weighted_trees.close();
        reweighted_input_trees.close();
        source_trees_for_this_search = "z_weighted_input_trees_reweighted";    
      }else {
        source_trees_for_this_search = argv[1];
      }




      int the_best_wrf_distance_of_current_hill= INT_MAX;
      string the_best_supertree_of_current_hill;

      //search until reaching a local optimum
      int iteration = 0;
      while (true) {
          iteration ++;

          //cout << "==========================" << init_supertree << endl;
          Node* myTree= build_tree(init_supertree);
          adjustTree(myTree);

          int total_nodes=0;
          total_number_of_nodes(myTree, total_nodes);


          //write init_supertree into "z_spr_neighbours"
          const char* neighbors_file = "z_spr_neighbours";
          writeToFile(init_supertree, neighbors_file);


          //writes all valid spr_neighbors into the file "z_spr_neighbours"
          //int number_of_neighbors= produce_all_spr_neighbors(myTree, total_nodes);
          int number_of_neighbors= produce_n_percent_of_spr_neighbors(myTree, total_nodes, neighborhood_percentage);
          //cout << "The number of neighbors in iteration " << iteration << " is: " << number_of_neighbors << endl;


          int best_distance_of_current_iter= INT_MAX;
          string best_supertree_of_current_iter;


          string suptree;

          ifstream supertrees;
          supertrees.open ("z_spr_neighbours");

          //This while goes through all the ST's in the
          // "z_spr_neighbours" file.
          while (getline(supertrees, suptree)) {

            //replace "suptree" with the ST at the beginning of the "z_wrf_input"
                //This is done in 3 steps:
                // 1- creating a new file, "z_wrf_input", in each iteration of the while-loop
                // 2- adding the "suptree" string variable at the beginning of "z_wrf_input"
                // 3- concatenating the source_trees' file to "z_wrf_input"

            std::ofstream wrf_file;
            wrf_file.open ("z_wrf_input");
            string current_supertree= suptree+"\n";
            wrf_file << current_supertree;

            ifstream source_trees(source_trees_for_this_search);
            wrf_file << source_trees.rdbuf();
            source_trees.close();

            wrf_file.close();


            //the following line computes the wrf dist. of ST to source trees
            int current_supertree_wrf_distance;
            if(ratchet) {
              current_supertree_wrf_distance = calculate_total_wrf("z_wrf_input");
              cout << "weightedrf: " << current_supertree_wrf_distance << endl;
            } else {
              current_supertree_wrf_distance = calculate_total_rf("z_wrf_input");
              cout << "rf: " <<current_supertree_wrf_distance << endl;
            }    
            
            ///cout << suptree << "  with score: " << current_supertree_wrf_distance << endl;
            if(current_supertree_wrf_distance < best_distance_of_current_iter) {
                best_distance_of_current_iter= current_supertree_wrf_distance;
                best_supertree_of_current_iter= suptree;
            }
          }//end of one local search
          supertrees.close();



      cout << "****************************************************" << endl;
      cout << "The bets ST found at the end of " << iteration << "th iteration and among "<<
          number_of_neighbors << " spr-neighbors is" <<endl;
      cout << best_supertree_of_current_iter << endl;
      cout << "And its WRF distance is: " << best_distance_of_current_iter << endl;
      cout << "****************************************************" << endl;

      //init_supertree = best_supertree_of_current_iter;  //it's WRONG to be here.
      /****************************
      NOTE:
      the above line, should be moved to inside of next if-statement CUZ when we reach a local opt,
      "best_supertree_of_current_iter" is not better than "init_supertree" at the beginning of this iteration,
      and we do not want to change "init_supertree". Note this does NOT matter for basi hill climbing, qsdt. 
      However, in ratchet, it DOES matter since the best supertree tree of each hill, will be fed as initial
      supertree of next hill climbing phase.
      *****************************/

      //Is it local optimum? 
      //update "the_best_wrf_distance_of_current_hill" if needed
      if(best_distance_of_current_iter < the_best_wrf_distance_of_current_hill) {

          init_supertree = best_supertree_of_current_iter;

          the_best_wrf_distance_of_current_hill= best_distance_of_current_iter;
          the_best_supertree_of_current_hill= best_supertree_of_current_iter;
      }else { // local optimum

          cout << "----------------------------------------------------------------------------------" << endl;
          
          if(ratchet) {
              cout << "ratchet local opt (re-weighted) reached." << endl;
          }else {
              cout << "regular local opt (original weights) reached. ###" << endl;
          }

          cout << "We have reached a local optimum which has better \(or equal\) WRF distance than all its spr-neighbors\n";
          cout << "The best SuperTree found after " << iteration << " number of iterations is: " << endl;
          cout << "And its WRF distance is " << the_best_wrf_distance_of_current_hill << endl;
          cout << "----------------------------------------------------------------------------------" << endl;
          cout << the_best_supertree_of_current_hill << endl;
          cout << "----------------------------------------------------------------------------------" << endl;


          
          if(!ratchet) {  //we are at the end of one ratchet iteration

              cout << "=======================================end of " << ratchet_counter << "-th ratchet iter=========================================" << endl;
              cout << "========================================================================================================" << endl;

              if(the_best_wrf_distance_of_current_hill < the_best_wrf_distance_seen){ //keep track of best supertree seen so far
                  the_best_wrf_distance_seen= the_best_wrf_distance_of_current_hill;
                  the_best_supertree_seen= the_best_supertree_of_current_hill;
              }                   

              break;  //end of second phase of ONE ratchet iteration
                             
          }else { //now perform a regular branch swapping
              the_best_wrf_distance_of_current_hill= INT_MAX;   //this line so necessary!!
                                                              //Because we are about to start a completely new hill climbing with new objective function.
                                                              //without this line, after 2nd ratchet iter, no improvement will be made since the weighted 
                                                              //distance is always larger than un-weighted and the init ST from previous step is local opt
                                                              //already.

              ratchet= false;
              iteration= 0;
          }

      }

  }//end of one branch swapping search loop


  }//end of ratchet search loop


  finish_time=clock();
  float diff ((float)finish_time - (float)start_time);
  float seconds= diff / CLOCKS_PER_SEC;
  //cout << "The #of clock ticks of this iteration: " << diff << endl;


  cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
  cout << "The best SuperTree found after " << number_of_ratchet_iterations << " number of ratchet iterations is: " << endl;    
  cout << "And its RF distance is " << the_best_wrf_distance_seen << endl;
  cout << "the running time is: " << seconds << " sec." << endl;        
  cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
  cout << the_best_supertree_seen << endl;
  cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
  cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$+" << endl;

  return 0;
}




//the following function was main() in original hashRF
//i create fake argv for it in main
double calculate_wrf_btwn_ST_n_source_tree(int argc, char** argv)
{
  double wrf_distance= 0.0;

  string outfilename;
  string infilename;
  bool bUbid = false; // for counting the number of unique bipartitions

  // TCLAP
  try {

    // Define the command line object.
    string 	helpMsg  = "HashRF\n";

    helpMsg += "Input file: \n";
    helpMsg += "   The current version of HashRF only supports the Newick format.\n";

    helpMsg += "Example of Newick tree: \n";
    helpMsg += "   (('Chimp':0.052625,'Human':0.042375):0.007875,'Gorilla':0.060125,\n";
    helpMsg += "   ('Gibbon':0.124833,'Orangutan':0.0971667):0.038875);\n";
    helpMsg += "   ('Chimp':0.052625,('Human':0.042375,'Gorilla':0.060125):0.007875,\n";
    helpMsg += "   ('Gibbon':0.124833,'Orangutan':0.0971667):0.038875);\n";

    helpMsg += "Print out mode: (defualt = matrix).\n";
    helpMsg += "   -p no, no output (default).\n";
    helpMsg += "   -p list, print RF distance in list format.\n";
    helpMsg += "   -p rate, print RF distance rate in list format.\n";
    helpMsg += "   -p matrix, print reuslting distance in matrix format.\n";

    helpMsg += "File option: \n";
    helpMsg += "   -o <export-file-name>, save the distance result in a file.\n";

    helpMsg += "Weighted RF distance: to select RF distance mode between weighted and unweighted (defualt = unweighted).\n";
    helpMsg += "   -w, compute weighted RF distance.\n";

    helpMsg += "Specify c value: \n";   
    helpMsg += "   -c <rate>, specify c value (default: 1000) \n";
    			
    helpMsg += "Examples: \n";
    helpMsg += "  hashf foo.tre 1000\n";
    helpMsg += "  hashf bar.tre 1000 -w\n";
    helpMsg += "  hashf bar.tre 1000 -w -p matrix\n";
    helpMsg += "  hashf bar.tre 1000 -w -p list\n";
    helpMsg += "  hashf bar.tre 1000 -w -p list -o output.dat\n";

    TCLAP::CmdLine cmd(helpMsg, ' ', "6.0.0");

    TCLAP::UnlabeledValueArg<string>  fnameArg( "name", "file name", true, "intree", "Input tree file name"  );
    cmd.add( fnameArg );

    TCLAP::UnlabeledValueArg<int>  numtreeArg( "numtree", "number of trees", true, 2, "Number of trees"  );
    cmd.add( numtreeArg );

    TCLAP::SwitchArg weightedSwitch("w", "weighted", "Compute weighted RF distance", false);
    cmd.add( weightedSwitch );

    TCLAP::ValueArg<string> printArg("p", "printoptions", "print options", false, "matrix", "Print options");
    cmd.add( printArg );

    TCLAP::ValueArg<unsigned int> cArg("c", "cvalue", "c value", false, 1000, "c value");
    cmd.add( cArg );
    
    TCLAP::SwitchArg ubidSwitch("u", "uBID", "unique BID count", false);
    cmd.add( ubidSwitch );
   
    TCLAP::ValueArg<string> outfileArg("o", "outfile", "Output file name", false, "", "Output file name");
    cmd.add( outfileArg );
    
 		// debug
		TCLAP::ValueArg<int> seedArg("s", "seedvalue", "user specified seed value", false, 1000, "user specified seed value");
    cmd.add( seedArg );
				
    cmd.parse( argc, argv );

    NUM_TREES = numtreeArg.getValue();

    if (NUM_TREES == 0) {
      string strFileLine;
      unsigned long ulLineCount;
      ulLineCount = 0;

      ifstream infile(argv[1]);
      if (infile) {
        while (getline(infile, strFileLine)) {
          ulLineCount++;
        }
      }
      //cout << "*** Number of trees in the input file: " << ulLineCount << endl;
      NUM_TREES = ulLineCount;

      infile.close();
    }

    if (NUM_TREES < 2) {cerr << "Fatal error: at least two trees expected.\n"; exit(2);}

    //if (weightedSwitch.getValue())
      //WEIGHTED = true;
    ///*************************************************************************
    //FOR SOME REASON above condition does NOT work correctly, and I replace it with this:
    if(argc == 4)  
      WEIGHTED = true;
    else if (argc == 3)
      WEIGHTED = false;

    if (printArg.getValue() != "matrix") {
      string printOption = printArg.getValue();
      if (printOption == "no")
        PRINT_OPTIONS = 0;
      if (printOption == "list")
        PRINT_OPTIONS = 1;
      if (printOption == "rate")
        PRINT_OPTIONS = 2;
      if (printOption == "cmb")
        PRINT_OPTIONS = 4;
    }

    if (cArg.getValue())
      C = cArg.getValue();
      
    if (seedArg.getValue())
      NEWSEED = seedArg.getValue();      
    
    if (ubidSwitch.getValue())
      bUbid = ubidSwitch.getValue();  
 
    outfilename = outfileArg.getValue();

  } catch (TCLAP::ArgException &e) { // catch any exceptions
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
  }


  /*****************************************************/
//  cout << "*** Reading a tree file and parsing the tree for taxon label collection ***\n";
  /*****************************************************/

  NEWICKTREE *newickTree;
  int err;
  FILE *fp;
  fp = fopen(argv[1], "r");
  if (!fp) { cout << "Fatal error: file open error\n"; exit(0); }


  newickTree = loadnewicktree2(fp, &err);
  if (!newickTree) {
    switch (err) {
    case -1:
      printf("Out of memory\n");
      break;
    case -2:
      printf("parse error\n");
      break;
    case -3:
      printf("Can't load file\n");
      break;
    default:
      printf("Error %d\n", err);
    }
  }

  /*****************************************************/
  //cout << "\n*** Collecting the taxon labels ***\n";
  /*****************************************************/
  LabelMap lm;

  try	{
    GetTaxaLabels(newickTree->root, lm);
  }
  catch (LabelMap::AlreadyPushedEx ex) { 
    cerr << "Fatal error: The label '" << ex.label << "' appeard twice in " << endl;
    exit(2);
  }
  NUM_TAXA = lm.size();
  //cout << "    Number of taxa = " << lm.size() << endl;
  killnewicktree(newickTree);
  fclose(fp);


  /*****************************************************/
  //cout << "\n*** Reading tree file and collecting bipartitions ***\n";
  /*****************************************************/
  HashRFMap vvec_hashrf; // Class HashRFMap

  ////////////////////////////
  // Init hash function class
  ////////////////////////////
  unsigned long long M1=0;
  unsigned long long M2=0;

	if (NEWSEED != 1000)
  	vvec_hashrf.uhashfunc_init(NUM_TREES, NUM_TAXA, C, NEWSEED);
  else
	 	vvec_hashrf.uhashfunc_init(NUM_TREES, NUM_TAXA, C);
 
  M1 = vvec_hashrf._HF.getM1();
  M2 = vvec_hashrf._HF.getM2();
  vvec_hashrf._hashtab2.resize(M1);

  fp = fopen(argv[1], "r");
  if (!fp) {cout << "Fatal error: file open error\n";  exit(0);}

  for (unsigned int treeIdx=0; treeIdx<NUM_TREES; ++treeIdx) {

    newickTree = loadnewicktree2(fp, &err);
    if (!newickTree) {
      switch (err) {
      case -1:
        printf("Out of memory\n");
        break;
      case -2:
        printf("parse error\n");
        break;
      case -3:
        printf("Can't load file\n");
        break;
      default:
        printf("Error %d\n", err);
      }
    }
    else {
      unsigned int numBitstr=0;

      dfs_compute_hash(newickTree->root, lm, vvec_hashrf, treeIdx, numBitstr, M1, M2);

      killnewicktree(newickTree);
    }
  }

  //cout << "    Number of trees = " << NUM_TREES << endl;
  fclose(fp);


  /*****************************************************/
  //cout << "\n*** Compute distance ***\n";
  /*****************************************************/
//  vector< vector<unsigned short> > SIM(NUM_TREES, vector<unsigned short>(NUM_TREES, 0)); // similarity matrix
  typedef vector< vector<float> > FLOAT_MATRIX_T;
  typedef vector< vector<unsigned short> > SHORT_MATRIX_T;
  SHORT_MATRIX_T SHORTSIM;
  FLOAT_MATRIX_T FLOATSIM; 

  //---------------------------------------- UNWEIGHTED --------------------------------------------
  if (!WEIGHTED) //  for unweighted RF
    SHORTSIM = SHORT_MATRIX_T (NUM_TREES, vector<unsigned short>(NUM_TREES,0));
  else // for weighted RF
    FLOATSIM = FLOAT_MATRIX_T (NUM_TREES, vector<float>(NUM_TREES,0.0));
       
	unsigned long uBID = 0;
	
	if (!WEIGHTED) {    // unweighted
    for (unsigned int hti=0; hti<vvec_hashrf._hashtab2.size(); ++hti) {
      unsigned int sizeVec = vvec_hashrf._hashtab2[hti].size();
			
			if (sizeVec) {

      	uBID += sizeVec;
      	
      	if (!bUbid) {
  	      for (unsigned int i=0; i<sizeVec; ++i) {
  	        unsigned int sizeTreeIdx = vvec_hashrf._hashtab2[hti][i]._vec_treeidx.size();
   
  	        if (sizeTreeIdx > 1) {
  	          for (unsigned int j=0; j<sizeTreeIdx; ++j) {
  	            for (unsigned int k=0; k<sizeTreeIdx; ++k) {
  	              if (j == k) continue;
  	              else {
  	                SHORTSIM[vvec_hashrf._hashtab2[hti][i]._vec_treeidx[j]][vvec_hashrf._hashtab2[hti][i]._vec_treeidx[k]] += 1;
  	              }
  	            }
  	          }
  	        }
  	      }
	      } //if
	    }
    }
  } 
  //---------------------------------------- WEIGHTED --------------------------------------------  
  else {
    vvec_hashrf._hashtab.resize(vvec_hashrf._hashtab2.size());
    for (unsigned int hti=0; hti<vvec_hashrf._hashtab2.size(); ++hti) {
      unsigned int sizeLinkedList = vvec_hashrf._hashtab2[hti].size();
      if (sizeLinkedList > 0) {
        for (unsigned int i1=0; i1<sizeLinkedList; ++i1) {
          unsigned int bidi = vvec_hashrf._hashtab2[hti][i1]._vec_treeidx.size();
          for (unsigned int i2=0; i2<bidi; ++i2) {
            BUCKET_STRUCT_T bk;
            bk._hv2 = vvec_hashrf._hashtab2[hti][i1]._hv2;
            bk._t_i = vvec_hashrf._hashtab2[hti][i1]._vec_treeidx[i2];
            bk._dist = vvec_hashrf._hashtab2[hti][i1]._vec_dist[i2];
            vvec_hashrf._hashtab[hti].push_back(bk);
          }
        }
      }
    }
    vvec_hashrf._hashtab2.clear();

    for (unsigned int hti=0; hti<vvec_hashrf._hashtab.size(); ++hti) {

      unsigned int sizeLinkedList = vvec_hashrf._hashtab[hti].size();

      if (sizeLinkedList > 1) {
        vector<unsigned long> vec_hv2;
        vector<unsigned long>::iterator itr_vec_hv2;

        ////////////////////////////////////////////////
        // Collect unique hv2 values in the linked list
        ////////////////////////////////////////////////
        for (unsigned int i=0; i<sizeLinkedList; ++i) {
          unsigned long hv2 = vvec_hashrf._hashtab[hti][i]._hv2;
          if (vec_hv2.empty())
            vec_hv2.push_back(hv2);
          else {
            itr_vec_hv2 = find(vec_hv2.begin(), vec_hv2.end(), hv2);
            if (itr_vec_hv2 == vec_hv2.end())
              vec_hv2.push_back(hv2);
          }
        }

        // distance array
        vector< vector<float> > vvec_dist(vec_hv2.size(), vector<float>(NUM_TREES, 0));

        /////////////////////////////////////////////////////////////
        // SET THE distance array with distance at proper tree index
        /////////////////////////////////////////////////////////////
        for (unsigned int i=0; i<sizeLinkedList; ++i) {
          for (unsigned int j=0; j<vec_hv2.size(); ++j) {
            if (vvec_hashrf._hashtab[hti][i]._hv2 == vec_hv2[j])
              vvec_dist[j][vvec_hashrf._hashtab[hti][i]._t_i] = vvec_hashrf._hashtab[hti][i]._dist;
          }
        }

        /////////////////////////////////////
        // UPDATE FLOATSIM MATIRX USING vvec_dist
        /////////////////////////////////////
        for (unsigned int i=0; i<vvec_dist.size(); ++i) {

          for (unsigned int j=0; j<vvec_dist[i].size(); ++j) {
            for (unsigned int k=0; k<vvec_dist[i].size(); ++k) {
              if (j == k) continue;
              else
                FLOATSIM[j][k] += abs(vvec_dist[i][j] - vvec_dist[i][k]);
            }
          }
        }

        vec_hv2.clear();
        vvec_dist.clear();
      }
      else if (sizeLinkedList == 1) {

        /////////////////////////////////////////////////////
        // PROPAGATE the dist value TO OTHER TREES' distance
        /////////////////////////////////////////////////////
        for (unsigned int i=0; i<NUM_TREES; ++i) {

          if (i == vvec_hashrf._hashtab[hti][0]._t_i)
            continue;
          else {
            FLOATSIM[i][vvec_hashrf._hashtab[hti][0]._t_i] += vvec_hashrf._hashtab[hti][0]._dist;
            FLOATSIM[vvec_hashrf._hashtab[hti][0]._t_i][i] += vvec_hashrf._hashtab[hti][0]._dist;
          }
        }
      }
    }
  }

  //cout << "    # of unique BIDs = " << uBID << endl;

	if (!WEIGHTED)
  	wrf_distance = print_rf_short_matrix(SHORTSIM, PRINT_OPTIONS, outfilename);
  else 
  	wrf_distance = print_rf_float_matrix(FLOATSIM, PRINT_OPTIONS, outfilename);


  /*****************************************************/
//  cout << "\n*** Print statistics ***\n";
  /*****************************************************/
  // CPU time comsumed
  struct rusage a;
  if (getrusage(RUSAGE_SELF,&a) == -1) { cerr << "Fatal error: getrusage failed.\n";  exit(2); }
  //cout << "\n    Total CPU time: " << a.ru_utime.tv_sec+a.ru_stime.tv_sec << " sec and ";
  //cout << a.ru_utime.tv_usec+a.ru_stime.tv_usec << " usec.\n";
 

  /*****************************************************/
//  cout << "\n*** Clear allocated memory***\n";
  /*****************************************************/
  SHORTSIM.clear();
  FLOATSIM.clear();
  vvec_hashrf.hashrfmap_clear();


  return wrf_distance;
}


// eof








//////////////////////////////////////REZA////////////////////////////////////>>>>>>>>>>>>>
//////////////////////////////////////////////////////////////////////////////>>>>>>>>>>>>>




//this function changes "percentage_to_be_reweighted"%, between 0 and 100, of weights to "new_weight". This is done by
//iterating through each weight usin regex, and based on some probability change it to "new_weight"
//NOTE there are two assumptions: "tree" is weighted already, AND "new_weight" should be only one digit
//Don't forget "tree" to be the original one with all weights = 1 (z_weighted_input_trees_original file).
//If you pass a tree which has been already re-weighted, then it is messed up, right? :D  
string change_branch_lengths(const string &tree, int percentage_to_be_reweighted, int new_weight){
    string re_weighted_tree = tree;
    char weight_char = '0'+new_weight;    //"i +'0'" is for conversion of one digit to char
                                          //NOTE: weiht should be one digit!!! 

    //find the location of weights, i.e. "1"s in "..):1", in newick format
    char befor_weight= ':';
    vector<int> weight_locations;
    for(int i =0; i < tree.size(); i++)
        if(tree[i] == befor_weight){
            weight_locations.push_back(i+1);
        }    

    //randomly select "percentage_to_be_reweighted"% of weights. We will reweight them to "new_weight"
    //here is how I do this:
    //1- create an array "arr" in which element i contains i
    //2- shuffle elements
    //3- give the first m elements weight "new_weight", where m = "percentage_to_be_reweighted"

    //step 1
    int num_internal_edges = weight_locations.size();
    int arr[weight_locations.size()];
    for(int i = 0; i < num_internal_edges; ++i){
        arr[i] = i;
    }
    
    //step 2    
    random_shuffle(arr, arr + num_internal_edges);

    //find m in step 3
    int num_of_edges_to_be_reweighted= (num_internal_edges*percentage_to_be_reweighted)/100;
    
    //step 3
    for(int i = 0; i < num_of_edges_to_be_reweighted; i++){
        re_weighted_tree[weight_locations[arr[i]]]= weight_char;
    }

    return re_weighted_tree;

}



//writes "s\n" into file "file_name"
void write_line_to_File(string s,const char* file_name) {

    ofstream myFile;
    myFile.open(file_name, std::ios_base::app);
    myFile << s + ";\n";
    myFile.close();

}


//given ST and one source tree in input_file, writes "restricted" ST and source tree in file output_file
//assumes ST and source tree are given in file "z_st_and_source_tree" and writes 'restricted' ST and source tree on "z_pruned_st_and_the_source_tree"
//it is too hard wired, I know, but should be fine I need it just to work now :))
void restrict_st_to_source_tree() {
/*
For now I am using my python script and terminall commands using system calls, but in the future 
you can implemet it like:
1- read first and second line of the input tree file (i.e. ST and one of source trees)
2- call your function to find non-shared taxa
3- traverse the ST tree data structure, you can use GetTaxaLabels(), to remove uncommon taxa AND removing nodes with 1 child
4- copy this new (restricted) ST and source trees into a file which will be fed to the nexat call of loadnewicktree2()
*/


	//the python program, "prune_tree.py", takes file "z_st_and_source_tree" containing ST and one source tree,
	// and prunes ST to taxa set of source tree, AND adds branch lengths 1 to all edges in newick format, AND
	//writes them into file "z_pruned_st_and_the_source_tree".
	std::string command1 = "python prune_tree.py z_st_and_source_tree";
	system(command1.c_str());
	//there are some unnecessary branch lengths which should be removed:
	//1- the branch lengths for interlan edges are like "..(t1,t2)1:1". Idk what is the first "1" for. I remove it.
	string command2 = "perl -pe 's/(?<=\\))1//g' z_pruned_st_and_the_source_tree > z_temp; cat z_temp> z_pruned_st_and_the_source_tree; rm z_temp";
	system(command2.c_str());
	//2- having branch length for leaves in this context, i.e. weighted RF distance is unnecessary.
	//thus, anything in newick format which is of the form "(some_taxon:1,...." should be removed. We keep only
	//branch lengths after afer closing parenthesis.
	string command3 = "perl -pe 's/(?<=[a-zA-Z0-9]):1//g' z_pruned_st_and_the_source_tree > z_temp; cat z_temp> z_pruned_st_and_the_source_tree; rm z_temp";
  system(command3.c_str());
}


//assumes "input_file" contains ST (FIRST line) and all source trees
double calculate_total_wrf(const char* input_file){
  
  double total_wrf_dist = 0.0;

  //fake argv and argc to be passed to calculate_wrf_btwn_ST_n_source_tree()
  //note here is how hashrf is called when we have exactly two tree: "hashrf z_pruned_st_and_the_source_tree 2 -w"
  //Thus:
  int fake_argc = 4;
  char* fake_argv[4];
  fake_argv[0] = "hashrf";
  fake_argv[1] = "z_pruned_st_and_the_source_tree";
  fake_argv[2] = "2";
  fake_argv[3] = "-w";

  string supertree;  
  string tree;
  int counter = 1;
  ifstream input;
  input.open (input_file);
  while (getline(input, tree)) {
    if(counter == 1) {  //ST is first tree
      supertree = tree+"\n";   
      counter++;
    }else {
      string current_tree= tree+"\n";
      ofstream temp_file;
      temp_file.open ("z_st_and_source_tree");
      temp_file << supertree;
      temp_file << current_tree;
      temp_file.close();
      
      //assumes ST and source tree are given in file "z_st_and_source_tree" and writes
      //'restricted' ST and source tree on "z_pruned_st_and_the_source_tree"
      restrict_st_to_source_tree();

      double dist= calculate_wrf_btwn_ST_n_source_tree(fake_argc, fake_argv);
      //cout << dist << endl;
      total_wrf_dist += dist;

    }  
  }
  input.close();

  return total_wrf_dist;

}

//assumes "input_file" contains ST (FIRST line) and all source trees
//SAME FUNCTION AS calculate_total_wrf(), BAD SMELL :||
double calculate_total_rf(const char* input_file){
  
  double total_rf_dist = 0.0;

  //fake argv and argc to be passed to calculate_wrf_btwn_ST_n_source_tree()
  //note here is how hashrf is called when we have exactly two tree: "hashrf z_pruned_st_and_the_source_tree 2 -w"
  //Thus:
  int fake_argc = 3;
  char* fake_argv[3];
  fake_argv[0] = "hashrf";
  fake_argv[1] = "z_pruned_st_and_the_source_tree";
  fake_argv[2] = "2";

  string supertree;  
  string tree;
  int counter = 1;
  ifstream input;
  input.open (input_file);
  while (getline(input, tree)) {
    if(counter == 1) {  //ST is first tree
      supertree = tree+"\n";   
      counter++;
    }else {
      string current_tree= tree+"\n";
      ofstream temp_file;
      temp_file.open ("z_st_and_source_tree");
      temp_file << supertree;
      temp_file << current_tree;
      temp_file.close();
      
      //assumes ST and source tree are given in file "z_st_and_source_tree" and writes
      //'restricted' ST and source tree on "z_pruned_st_and_the_source_tree"
      restrict_st_to_source_tree();

      double dist= calculate_wrf_btwn_ST_n_source_tree(fake_argc, fake_argv);
      //cout << dist << endl;
      total_rf_dist += dist;

    }  
  }
  input.close();

  return total_rf_dist;

}








//*******************added for SPR neighborhood*********************
//******************************************************************
//applies all the possible spr's on "tree", and writes them into the file "z_spr_neighbours"
int produce_all_spr_neighbors(Node* myTree, int number_of_taxa) {

    int number_of_neighbors= 0;


    //The following line deletes the contents of z_spr_neighbours
    const char* neighbors_file = "z_spr_neighbours";
    std::ofstream ofile(neighbors_file, ios_base::trunc);

    for(int i=1; i<number_of_taxa; i++) {

        Node* spr_on= myTree->find_by_prenum(i);
        int which_sibling=0;
        /*
        cout << "************************************************************" << endl;
        cout << "spr_on's pre-order number = i = " << spr_on->get_preorder_number() << endl;
        cout << "new_sibling's pre-order number = j = " << spr_on->get_preorder_number() << endl;
        cout << "original  tree: " << myTree->str_subtree() << endl;
        cout << "************************************************************" << endl;
        */

        for(int j=1; j<number_of_taxa; j++) {

            if (j != i) {

                Node* new_sibling= myTree->find_by_prenum(j);
                bool good_spr= true;
                //cout<< "i = " <<  i << ", j = " << j << endl;



                //bad spr: check whether new_sibling is parent
                Node* parent= spr_on->get_p();
                if (parent->get_preorder_number() == j) {
                    //cout << "-----BAD SPR WAS IGNORED----- \n" << endl;
                    continue;
                }


                //bad spr: check whether new_sibling is a descendant of spr_on
                if(! spr_on->is_leaf()) {
                    vector<Node *> descendants = spr_on->find_descendants();
                    vector<Node *>::iterator it;
                    for(it = descendants.begin(); it!= descendants.end(); ++it) {
                        if(new_sibling->get_preorder_number() == (*it)->get_preorder_number()) {
                            //cout << "-----BAD SPR WAS IGNORED----- \n" << endl;
                            good_spr= false;
                            continue;       //goes out of the inner loop which is 4 lines above
                        }
                    }
                }


                

                //bad spr: if new sibling is old sibling ignore it cuz this spr-move results to the original tree
                list<Node *>::iterator itr;
                for(itr = (spr_on->parent())->get_children().begin(); itr!= (spr_on->parent())->get_children().end(); ++itr) {
                    if((*itr)->get_preorder_number() == j) {
                        good_spr= false;
                        continue;
                    }
                }


                if(good_spr) {
                  Node* undo= spr_on->spr(new_sibling, which_sibling);
                  adjustTree(myTree);
                  number_of_neighbors ++;
                  string new_tree= myTree ->str_subtree();
                  writeToFile(new_tree, neighbors_file);
                  //cout<<  "tree after spr: " << new_tree << "\n" <<endl;

                  //restore the original tree
                  spr_on->spr(undo, which_sibling);
                  adjustTree(myTree);
              }

          }

        }

    }



//    cout << "**********************************************************" << endl;
//    cout << "Total Number Of Taxa In Trees: " << number_of_taxa << endl;
//    cout << "Current Node (Supertree) : " << myTree->str_subtree() << endl;
//    cout << "#of spr_neighbors of initial supertree = " << number_of_neighbors << endl;
//    cout << "Trees were written in z_spr_neighbours in newick format." << endl;
//    cout << "**********************************************************" << endl;
    return number_of_neighbors;

}



int produce_n_percent_of_spr_neighbors(Node* myTree, int number_of_taxa, int percentage) {

    int number_of_neighbors= 0;

    //this is the idea how I produce a subset of neighbors:
    //For producing all spr-neighbors we iterate over all nodes by their pre-order counts as "spr_on" node,
    //and then produce all possible valid neighbors. Now, I do the following to produce a subset of neghbors:
    //1- given #of taxa and percentage, calculate expected number of neighbors to be produced, M
    //2- produce M random numbers in range 1 and #of taxa (this is done the way I did it for ratchet_weights), call this set of nodes subset_of_nodes
    //3- for each node in subset_of_nodes, produce all possible and valid spr neighbors
    //NOTE that this approach may produce less OR more number of neghbors than |all_spr_neighbors|*percentage

    
    
    int arr[number_of_taxa];
    for(int i = 0; i < number_of_taxa; ++i)
        arr[i] = i;
    
    random_shuffle(arr, arr + number_of_taxa);
    
    int M = (number_of_taxa*percentage)/100;
    int subset_of_nodes[M];

    for(int i=0; i<M; i++) {
        subset_of_nodes[i]= arr[i];
    }



    //The following line deletes the contents of z_spr_neighbours
    const char* neighbors_file = "z_spr_neighbours";
    std::ofstream ofile(neighbors_file, ios_base::trunc);

    for(int counter=0; counter<M; counter++) {

        int i = subset_of_nodes[counter];
        if (i == 0) continue;   //SPR on root?? Nonsense.

        Node* spr_on= myTree->find_by_prenum(i);
        int which_sibling=0;
        /*
        cout << "************************************************************" << endl;
        cout << "spr_on's pre-order number = i = " << spr_on->get_preorder_number() << endl;
        cout << "new_sibling's pre-order number = j = " << spr_on->get_preorder_number() << endl;
        cout << "original  tree: " << myTree->str_subtree() << endl;
        cout << "************************************************************" << endl;
        */

        for(int j=1; j<number_of_taxa; j++) {

            if (j != i) {

                Node* new_sibling= myTree->find_by_prenum(j);
                bool good_spr= true;
                //cout<< "i = " <<  i << ", j = " << j << endl;



                //bad spr: check whether new_sibling is parent
                Node* parent= spr_on->get_p();
                if (parent->get_preorder_number() == j) {
                    //cout << "-----BAD SPR WAS IGNORED----- \n" << endl;
                    continue;
                }


                //bad spr: check whether new_sibling is a descendant of spr_on
                if(! spr_on->is_leaf()) {
                    vector<Node *> descendants = spr_on->find_descendants();
                    vector<Node *>::iterator it;
                    for(it = descendants.begin(); it!= descendants.end(); ++it) {
                        if(new_sibling->get_preorder_number() == (*it)->get_preorder_number()) {
                            //cout << "-----BAD SPR WAS IGNORED----- \n" << endl;
                            good_spr= false;
                            continue;       //goes out of the inner loop which is 4 lines above
                        }
                    }
                }


                

                //bad spr: if new sibling is old sibling ignore it cuz this spr-move results to the original tree
                list<Node *>::iterator itr;
                for(itr = (spr_on->parent())->get_children().begin(); itr!= (spr_on->parent())->get_children().end(); ++itr) {
                    if((*itr)->get_preorder_number() == j) {
                        good_spr= false;
                        continue;
                    }
                }


                if(good_spr) {
                  Node* undo= spr_on->spr(new_sibling, which_sibling);
                  adjustTree(myTree);
                  number_of_neighbors ++;
                  string new_tree= myTree ->str_subtree();
                  writeToFile(new_tree, neighbors_file);
                  //cout<<  "tree after spr: " << new_tree << "\n" <<endl;

                  //restore the original tree
                  spr_on->spr(undo, which_sibling);
                  adjustTree(myTree);
              }

          }

        }

    }



//    cout << "**********************************************************" << endl;
//    cout << "Total Number Of Taxa In Trees: " << number_of_taxa << endl;
//    cout << "Current Node (Supertree) : " << myTree->str_subtree() << endl;
//    cout << "#of spr_neighbors of initial supertree = " << number_of_neighbors << endl;
//    cout << "Trees were written in z_spr_neighbours in newick format." << endl;
//    cout << "**********************************************************" << endl;
    return number_of_neighbors;

}





//writes "s\n" into file "file_name"
void writeToFile(string s,const char* file_name) {

    ofstream myFile;
    myFile.open(file_name, std::ios_base::app);
    myFile << s + ";\n";
    myFile.close();

}


void adjustTree(Node* myTree) {

    myTree->set_depth(0);
    myTree->fix_depths();
    myTree->preorder_number();
    myTree->edge_preorder_interval();

}




//returns the total number of nodes in the clade of "node", i.e. (#of taxa)+(#of internal nodes)
//This method is a modification of the set_preorder_number() method in node.h
void total_number_of_nodes(Node* node, int& total_nodes) {
    total_nodes++;
    list<Node *>::iterator c;
    list<Node *> children= node->get_children();
    for(c = children.begin(); c != children.end(); c++) {
        total_number_of_nodes(*c, total_nodes);
    }
}



