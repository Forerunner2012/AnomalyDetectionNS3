/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */

#include "ns3/core-module.h"
#include "ns3/ad-detector-simulator-helper.h"

#include "ns3/mobility-module.h"
#include "ns3/ns2-mobility-helper.h"

using namespace ns3;
using namespace std;

NS_LOG_COMPONENT_DEFINE ("ReadCoordinatesToNs2MobilityHelper");
//export NS_LOG=MobilitySimulator=level_all

// ---------- End of Namespaces ----------------------------------------------

// ---------- Prototypes -----------------------------------------------------
static void CourseChange (std::ostream *os, std::string foo, Ptr<const MobilityModel> mobility);
// ---------- End of Prototypes ----------------------------------------------

int 
main (int argc, char *argv[])
{
    NS_LOG_UNCOND ("Mobility Simulator for Training");
    LogComponentEnable ("ReadCoordinatesToNs2MobilityHelper", LOG_LEVEL_ALL);
    // ---------- Simulation Variables -----------------------------------------
    // Change the variables and file names only in this block!
    //NOTE Put variable that are needed for the simulation
    std::string traceFile;
    std::string logFile;
    int    nodeNum;
    double duration;
    
    // Enable logging from the ns2 helper
    //LogComponentEnable ("MobilitySimulator",LOG_LEVEL_DEBUG);
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
        "./waf --visualize --run \"read-coordinatesToNs2MobilityHelper"
        " --traceFile=src/ad-detector-simulator/data/A6-d07-h08_ns2formated.log"
        " --nodeNum=2082 --duration=1800.0 --logFile=ns2-mob.log\" \n\n"
        "NOTE: ns2-traces-file could be an absolute or relative path. You could use the file default.ns_movements\n"
        "      included in the same directory of this example file.\n\n"
        "NOTE 2: Number of nodes present in the trace file must match with the command line argument and must\n"
        "        be a positive number. Note that you must know it before to be able to load it.\n\n"
        "NOTE 3: Duration must be a positive number. Note that you must know it before to be able to load it.\n\n";
        return 0;
    }
    // ---------- End of Parameter from Command line ---------------------------

    // ---------- Nodes Setup --------------------------------------------------
    //NOTE Setup the nodes (malicious nodes and non malicious ones)
    //NOTE think to not introduce malicious node during training
    NS_LOG_UNCOND ("Nodes Setup from traceFile: " << traceFile.c_str ());
    // Create Ns2MobilityHelper with the specified trace log file as parameter
    Ns2MobilityHelper ns2 = Ns2MobilityHelper (traceFile);

    // Open log file for output
    std::ofstream os;
    os.open (logFile.c_str ());
    
    // Create all nodes
    NS_LOG_DEBUG ("Creating nodes");
    NodeContainer stas;
    stas.Create (nodeNum);
    
    // Configure movements for each node, while reading trace file
    NS_LOG_DEBUG ("Installing nodes");
    ns2.Install ();
    // ---------- End of Nodes Setup -------------------------------------------
    
    // ---------- Network Setup ------------------------------------------------
    //NOTE Setup the network between nodes (Wave)
    NS_LOG_DEBUG ("Network Setup between nodes");
    
    // ---------- End of Network Setup -----------------------------------------

    // ---------- Allocate Node Positions --------------------------------------
    //TODO Setup the initialisation of the network
    // ---------- End of Allocate Node Positions -------------------------------

    // ---------- Use Mobibility Model for Node Positions ----------------------
    //TODO Affect Mobitlity to the node according to t
    // Configure callback for logging
    Config::Connect ("/NodeList/*/$ns3::MobilityModel/CourseChange", MakeBoundCallback (&CourseChange, &os));
    // ---------- End of Mobibility Model --------------------------------------


    // ---------- Simulation Monitoring ----------------------------------------
    //NS_LOG_INFO ("Configure Tracing.");
    //TODO: AsciiTraceHelper ascii
    
    NS_LOG_UNCOND ("Simulation is ready");
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

