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

#include <iostream>
#include <fstream>
#include <sstream>

#include "ns3/core-module.h"
#include "ns3/mobility-module.h"
#include "ns3/mobility-module.h"
#include "ns3/ns2-mobility-helper.h"

// ---------- Namespaces -----------------------------------------------------
using namespace ns3;
using namespace std;
NS_LOG_COMPONENT_DEFINE ("MobilitySimulator");
//export NS_LOG=MobilitySimulator=level_all
// ---------- End of Namespaces ----------------------------------------------

// ---------- Prototypes -----------------------------------------------------
static void CourseChange (std::ostream *os, std::string foo, Ptr<const MobilityModel> mobility);
// ---------- End of Prototypes ----------------------------------------------

int 
main (int argc, char *argv[])
{
    NS_LOG_UNCOND ("Mobility Simulator for Training");

    // ---------- Simulation Variables -----------------------------------------
    // Change the variables and file names only in this block!
    //NOTE Put variable that are needed for the simulation
    std::string traceFile;
    std::string logFile;
    int    nodeNum;
    double duration;
    
    // Enable logging from the ns2 helper
    LogComponentEnable ("Ns2MobilityHelper",LOG_LEVEL_DEBUG);
    // ---------- End of Simulation Variables ----------------------------------

    
    // ---------- Parameter from Command line ----------------------------------
    //NOTE Read and take into account commands lines from terminal
    // Parse command line attribute
    CommandLine cmd;
    cmd.AddValue ("traceFile", "Ns2 movement trace file", traceFile);
    cmd.AddValue ("nodeNum", "Number of nodes", nodeNum);
    cmd.AddValue ("duration", "Duration of Simulation", duration);
    cmd.AddValue ("logFile", "Log file", logFile);
    cmd.Parse (argc,argv);
    
    // Check command line arguments
    if (traceFile.empty () || nodeNum <= 0 || duration <= 0 || logFile.empty ())
    {
        std::cout << "Usage of " << argv[0] << " :\n\n"
        "./waf --run \"ns2-mobility-trace"
        " --traceFile=src/mobility/examples/default.ns_movements"
        " --nodeNum=2 --duration=100.0 --logFile=ns2-mob.log\" \n\n"
        "NOTE: ns2-traces-file could be an absolute or relative path. You could use the file default.ns_movements\n"
        "      included in the same directory of this example file.\n\n"
        "NOTE 2: Number of nodes present in the trace file must match with the command line argument and must\n"
        "        be a positive number. Note that you must know it before to be able to load it.\n\n"
        "NOTE 3: Duration must be a positive number. Note that you must know it before to be able to load it.\n\n";
        return 0;
    }
    // ---------- End of Parameter from Command line ---------------------------

    // ---------- Nodes Setup --------------------------------------------------
    //TODO Setup the nodes (malicious nodes and non malicious ones)
    //NOTE think to not introduce malicious node during training


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
    
    // Configure callback for logging
    Config::Connect ("/NodeList/*/$ns3::MobilityModel/CourseChange", MakeBoundCallback (&CourseChange, &os));
    Simulator::Stop (Seconds (duration));
    Simulator::Run ();
    
    
    Simulator::Destroy ();
    
    
    os.close (); // close log file
    // ---------- End of Simulation Monitoring ---------------------------------
    return 0;
}

// ---------- Function Definitions -------------------------------------------

//! Prints actual position and velocity when a course change event occurs
static void
CourseChange (std::ostream *os, std::string foo, Ptr<const MobilityModel> mobility) {
  Vector pos = mobility->GetPosition (); // Get position
  Vector vel = mobility->GetVelocity (); // Get velocity

  // Prints position and velocities
  *os << Simulator::Now () << " POS: x=" << pos.x << ", y=" << pos.y
      << ", z=" << pos.z << "; VEL:" << vel.x << ", y=" << vel.y
      << ", z=" << vel.z << std::endl;
}
// ---------- End of Tools functions Definitions -----------------------------

