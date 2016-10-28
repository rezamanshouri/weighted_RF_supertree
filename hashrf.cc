
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



static void
print_rf_short_matrix(
  vector< vector<unsigned short> > &SHORTSIM,
  unsigned options,
  string outfile)
{
  ofstream fout;
  if (outfile != "")
  {
    fout.open(outfile.c_str());
  }

  switch (options) {
  case 0:
    return;
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
    cout << "\nRobinson-Foulds distance (matrix format):\n";
    if (outfile == "") {
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
         wrf_dist =  float((SIM[0][1] + SIM[0][1])/4);
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

  //add weight 1 to all nternal edges: using perl replace all ")" with ""):1" except last one
  const char* input_file_weighted = "z_input_file_weighted";
  string argv1(argv[1]);
  string command1 = "perl -pe 's/\\)/):1/g' " + argv1 + " > z_temp; cat z_temp > " + input_file_weighted + "; rm z_temp";
  system(command1.c_str());
  string command2 = "perl -pe 's/:1;/;/g' z_input_file_weighted > z_temp; cat z_temp> z_input_file_weighted; rm z_temp";
  system(command2.c_str());
  

  double total_wrf_dist = calculate_total_wrf(argv[1]);
  cout << total_wrf_dist << endl;


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

    if (weightedSwitch.getValue())
      WEIGHTED = true;

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
  	print_rf_short_matrix(SHORTSIM, PRINT_OPTIONS, outfilename);
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
//NOTE: weiht should be one digit!!!
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
//too hard wired, I know. I need it just to work now :))
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





