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

#include <vector>
#include <string>

#include "ns3/core-module.h"
#include "ns3/mobility-module.h"

// ---------- Namespaces -----------------------------------------------------
namespace ns3{
    const bool DEBUG_MODE = true;
    const int N_DATA_COLUMN = 5;
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
// ---------- End of Prototypes ----------------------------------------------

int 
main (int argc, char *argv[])
{
    NS_LOG_UNCOND ("Mobility Simulator for Training");

    // ---------- Simulation Variables -----------------------------------------
    // Change the variables and file names only in this block!
    //NOTE Put variable that are needed for the simulation
    
    //std::string node_coordinates_file_name ("scratch/A6-d07-h08.log");
    std::string node_coordinates_file_name ("scratch/A6-d07-h08smaller.log");
    
    // ---------- End of Simulation Variables ----------------------------------

    
    // ---------- Parameter from Command line ----------------------------------
    //TODO Read and take into account commands lines from terminal
    CommandLine cmd;
    cmd.Parse (argc, argv);
    // ---------- End of Parameter from Command line ---------------------------

    // ---------- Nodes Setup --------------------------------------------------
    //TODO Setup the nodes (malicious nodes and non malicious ones)
    //NOTE think to not introduce malicious node during training
    NS_LOG_INFO ("Reading nodes coordonates from" << node_coordinates_file_name.c_str());
    vector<vector<double> > coord_array;
    coord_array = readCordinatesFile(node_coordinates_file_name);
    //printCoordinateArray(node_coordinates_file_name.c_str(), coord_array);    //DEBUG
    
    int n_nodes = max_element_column(coord_array,1); //NOTE identify position of data can be interesting to do with global constants)
    int n_time = max_element_column(coord_array,0);
    NS_LOG_INFO ("Number of nodes: " << n_nodes);
    NS_LOG_INFO ("Duration of the simulation: " << n_time <<"s");
    
    NS_LOG_INFO ("Create Nodes.");
    NodeContainer nodes;
    nodes.Create (n_nodes);
    // ---------- End of Nodes Setup -------------------------------------------
    
    // ---------- Network Setup ------------------------------------------------
    //TODO Setup the network (wave)
    // ---------- End of Network Setup -----------------------------------------

    // ---------- Allocate Node Positions --------------------------------------
    //TODO Setup the initialisation of the network
    // ---------- End of Allocate Node Positions -------------------------------

    // ---------- Use Mobibility Model for Node Positions ----------------------
    //TODO Affect Mobitlity to the node according to t
    
    //
    //TODO Introduce a function for each coursechange (see examples/mobitily/main-random-walk.cc with Config::Connect ("/NodeList/*/$ns3::MobilityModel/CourseChange",                    MakeCallback (&CourseChange));
    //
    
    // ---------- End of Mobibility Model --------------------------------------


    // ---------- Simulation Monitoring ----------------------------------------
    //NS_LOG_INFO ("Configure Tracing.");
    //TODO: AsciiTraceHelper ascii

    Simulator::Run ();

    // Collect data there

    Simulator::Destroy ();
    // ---------- End of Simulation Monitoring ---------------------------------
    return 0;
}

// ---------- Function Definitions -------------------------------------------
//! Read coordinate file (took and adapted from examples/matric-topology/matrix-topology.cc)
vector<vector<double> > readCordinatesFile (std::string node_coordinates_file_name)
{
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

        if (n != N_DATA_COLUMN) {
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
void printCoordinateArray (const char* description, vector<vector<double> > coord_array)
{
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
int max_element_column ( vector<vector<double> > coord_array, int column_number)
{
    if ((column_number > N_DATA_COLUMN) | (coord_array.at(0).size() > N_DATA_COLUMN))
        NS_LOG_ERROR ("ERROR: Number of columns provided is not equal in the expected range of value (0," << N_DATA_COLUMN << ")");
    
    int result = coord_array.at(0).at(column_number);
    for (size_t n = 1; n < coord_array.size (); n++)
    {
        if(result < coord_array.at(n).at(column_number))
            result = coord_array.at(n).at(column_number);
    }
    return result;
}
// ---------- End of Tools functions Definitions -----------------------------

