/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
To test:
./waf --visualize --run "ns2-mobility-trace \
--traceFile=scratch/A6-d07-h08_ns2formated.log \
--nodeNum=2083 \
--duration=1800.0 \
--logFile=ns2-mob.log"
*/

#include <vector>
#include <string>

#include "ns3/core-module.h"
#include "ns3/mobility-module.h"

// ---------- Namespaces -----------------------------------------------------
namespace ns3{
    const bool DEBUG_MODE = true;
    const int DATA_COLUMN_NUMBER = 5;
    const int DATA_COLUMN_TIME = 0;
    const int DATA_COLUMN_ID = 1;
    const int DATA_COLUMN_POSX = 2;
    const int DATA_COLUMN_POSY = 3;
    const int DATA_COLUMN_SPEED = 4;
}
using namespace ns3;
using namespace std;
NS_LOG_COMPONENT_DEFINE ("MobilitySimulator");
//export NS_LOG=MobilitySimulator=level_all
// ---------- End of Namespaces ----------------------------------------------

// ---------- Prototypes -----------------------------------------------------
vector<vector<double> > readCordinatesFile (std::string node_coordinates_file_name);
void printCoordinateArray (const char* description, vector<vector<double> > coord_array);
void printMatrix (const char* description, vector<vector<bool> > array);
int max_element_column (vector<vector<double> > coord_array, int column_number);
void convertCoordinateToFile (std::string convert_coordinates_file_name, vector<vector<double> > coord_array, int n_nodes);
// ---------- End of Prototypes ----------------------------------------------

int 
main (int argc, char *argv[])
{
    NS_LOG_UNCOND ("Tools for converting rawdata from Madrid to ns2 HelperMobility Model");


    
    // ---------- Simulation Variables -----------------------------------------
    // Change the variables and file names only in this block!
    //NOTE Put variable that are needed for the simulation
    
    //scratch/A6-d07-h08.log
    //std::string node_coordinates_file_name ("scratch/A6-d07-h08smaller.log");
    //convert_coordinates_file_name.insert(node_coordinates_file_name.size()-4,"_ns2formated");
    
    std::string node_coordinates_file_name;
    std::string convert_coordinates_file_name;
    // ---------- End of Simulation Variables ----------------------------------

    
    // ---------- Parameter from Command line ----------------------------------
    // Parse command line attribute
    CommandLine cmd;
    cmd.AddValue ("rawData", "Input data file for conversion", node_coordinates_file_name);
    cmd.AddValue ("convertData", "Data file converted to ns2", convert_coordinates_file_name);
    cmd.Parse (argc, argv);
    
    //Check command line arguments
    if (node_coordinates_file_name.empty () || convert_coordinates_file_name.empty ())
    {
        std::cout << "Usage of " << argv[0] << " :\n\n"
                "./waf --run \"scratch/tranform_rawdata_to_ns3mobilityhelper"
                " --rawData=scratch/A6-d07-h08.log"
                " --convertData=scratch/A6-d07-h08_ns2formated.log\" \n\n"
                "NOTE: ns2-traces-file could be an absolute or relative path. You could use the file default.ns_movements\n"
                "      included in the same directory of this example file.\n\n";
        return 0;
    }
    // ---------- End of Parameter from Command line ---------------------------

    
    // ---------- Conversion of Data  ------------------------------------------
    //TODO Setup the nodes (malicious nodes and non malicious ones)
    //NOTE think to not introduce malicious node during training
    NS_LOG_INFO ("Reading nodes coordonates from" << node_coordinates_file_name.c_str());
    vector<vector<double> > coord_array;
    coord_array = readCordinatesFile(node_coordinates_file_name);
    //printCoordinateArray(node_coordinates_file_name.c_str(), coord_array);    //DEBUG
    
    int n_nodes = max_element_column(coord_array,DATA_COLUMN_ID);
    int n_time = max_element_column(coord_array,DATA_COLUMN_TIME);
    cout << "Number of nodes detected: " << n_nodes << "(+1 if plan to use trace)" << endl;
    cout << "Duration of the simulation found: " << n_time <<"s" << endl;
        
    convertCoordinateToFile(convert_coordinates_file_name, coord_array, n_nodes);
    cout << "Converted data have been written in " << convert_coordinates_file_name << endl;
    // ---------- End of Conversion of Data  -----------------------------------
    
    
    // ---------- Simulation Monitoring ----------------------------------------
    Simulator::Run ();
    Simulator::Destroy ();
    // ---------- End of Simulation Monitoring ---------------------------------
    
    return 0;
}

// ---------- Function Definitions -------------------------------------------

//! Read coordinate file (took and adapted from examples/matric-topology/matrix-topology.cc)
vector<vector<double> > readCordinatesFile (std::string node_coordinates_file_name) {
    ifstream node_coordinates_file;
    node_coordinates_file.open (node_coordinates_file_name.c_str (), ios::in);
    
    if (node_coordinates_file.fail ())
    {
        NS_FATAL_ERROR ("File " << node_coordinates_file_name.c_str () << " not found");
    }
    
    vector<vector<double> > coord_array;
    int m = 0;
    while (!node_coordinates_file.eof ())
    {
        string line;
        getline (node_coordinates_file, line);

        if (line == "")
        {
            NS_LOG_WARN ("WARNING: Ignoring blank row: " << m);
            break;
        }

        istringstream iss (line);
        double coordinate;
        vector<double> row;
        int n = 0;
        while (iss >> coordinate)
        {
            row.push_back (coordinate);
            n++;
        }

        if (n != DATA_COLUMN_NUMBER) {
            NS_LOG_ERROR ("ERROR: Number of elements at line#" << m << " is "  << n << " which is not equal to 2 for node coordinates file");
            exit (-1);
        }
        else {
            coord_array.push_back (row);
        }
        m++;
    }
    node_coordinates_file.close ();
    return coord_array;
}

//! Print the coordinates loaded in vectors (took and adapted from examples/matric-topology/matrix-topology.cc)
void printCoordinateArray (const char* description, vector<vector<double> > coord_array) {
  cout << "**** Start " << description << "********" << endl;
  for (size_t m = 0; m < coord_array.size (); m++)
    {
      for (size_t n = 0; n < coord_array.at(m).size (); n++)
        {
          cout << coord_array.at(m).at(n) << ' ';
        }
      cout << endl;
    }
  cout << "**** End " << description << "********" << endl;
}

//! Look for the max in the given column
int max_element_column ( vector<vector<double> > coord_array, int column_number) {
    if ((column_number > DATA_COLUMN_NUMBER) | (coord_array.at(0).size() > DATA_COLUMN_NUMBER))
        NS_LOG_ERROR ("ERROR: Number of columns provided is not equal in the expected range of value (0," << DATA_COLUMN_NUMBER << ")");
    
    int result = coord_array.at(0).at(column_number);
    for (size_t n = 1; n < coord_array.size (); n++)
    {
        if(result < coord_array.at(n).at(column_number))
            result = coord_array.at(n).at(column_number);
    }
    return result;
}

//!
//! Write the coordinates converted to the file
void convertCoordinateToFile (std::string convert_coordinates_file_name, vector<vector<double> > coord_array, int n_nodes) {
    
    ofstream myfile;
    myfile.open (convert_coordinates_file_name.c_str());
    
    bool first_iteration[n_nodes];
    for(int i = 0; i < n_nodes; ++i)
		first_iteration[i] = false;
    
    //We need to treat each node one by one sequentially
    for (int k = 0; k < n_nodes; k++) {
        for (size_t n = 0; n < coord_array.size (); n++) {

            if (k == coord_array.at(n).at(DATA_COLUMN_ID)) {
                //Refer to the number we are dealing with (0, 1, 2 ...)

                //First time node appear? if y, need ti set pos x and y
                //$node_(0) set X_ 150.0
                //$node_(0) set Y_ 93.98597018956875
                if(first_iteration[k] == false) {
                    first_iteration[k] = true;
                    myfile << "$node_(" << k
                        << ") set X_ " << coord_array.at(n).at(DATA_COLUMN_POSX)
                        << endl;
                    myfile << "$node_(" << k
                        << ") set Y_ " << coord_array.at(n).at(DATA_COLUMN_POSY)
                        << endl;
                }

                //Basic line for describe the node
                //$ns_ at 0.0 "$node_(0) setdest 150.0 110.0 50.40378694202284"
                myfile << "$ns_ at " << coord_array.at(n).at(DATA_COLUMN_TIME) 
                    << " \"$node_(" << k << ") setdest "
                    << coord_array.at(n).at(DATA_COLUMN_POSX) << " "
                    << coord_array.at(n).at(DATA_COLUMN_POSY) << " "
                    << coord_array.at(n).at(DATA_COLUMN_SPEED) << "\""
                    << endl;
            }
        }     
    }
    cout << endl;
    myfile.close();
}
// ---------- End of Tools functions Definitions -----------------------------

