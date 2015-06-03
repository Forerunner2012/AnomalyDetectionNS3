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
/*! \file ads-simulator.cc
    \brief File containing the AD detector and simulation 
    This file contains both headers and declarations for classes + run the simulator in the main class (instead of classic couple .h and .cc because I need to understand the setup to do so in NS3 (act as module or similar with wscript).
*/

// ---------- Header & Includes ----------------------------------------------
#include "ns3/core-module.h"
#include <math.h>
#include "ns3/node.h"
/*
#include "ns3/log.h"
#include "ns3/synchronizer.h"
#include "ns3/object.h"
#include "ns3/object-factory.h"
#include "ns3/simulator.h"
#include "ns3/mobility-model.h"
#include "ns3/simulator-module.h"
#include "ns3/node-module.h"
#include "ns3/helper-module.h"
#include <math.h>
#include <iostream>
#include <cerrno>
#include <iostream>
#include <fstream>
#include <string>
*/

// ---------- Namespaces -----------------------------------------------------
using namespace ns3;
using namespace std;    //NOTE: a bit dangerous to do, caer about that later
// ---------- End of Namespaces ----------------------------------------------

// ---------- Prototypes -----------------------------------------------------
class Matrix;
class DataSource;
class Detector;
class EuclideanModel;
class DetectorEuclidean;
// ---------- End of Prototypes ----------------------------------------------

// ---------- Typedef --------------------------------------------------------
//! \enum Attacks types
//TODO Update for MC attacks
typedef enum {
    AT_DRAIN,
    AT_DRAIN_MIT,
    AT_GREYHOLE,
    AT_GREYHOLE_MIT,
    AT_FLOOD,
    AT_FLOOD_MIT,
    AT_MAXNUM} AttackType;

//! \typedef Structure to gather the result
typedef struct {
    bool decision;
    double believe;
    AttackType attack;
    } Result;
// ---------- End of Typedef -------------------------------------------------

// ---------- Classes --------------------------------------------------------
/*!
 * \class Matrix
 * \brief Tool class used for manipulate data from input data and features
 */
class Matrix {
private:
  int sizeI, sizeJ;
  double **pMatrix;

public:
    Matrix(int i, int j) {
        pMatrix = new double*[j];
        for(int it=0;it<j;++it) pMatrix[it] = new double[i];
        sizeI = i;
        sizeJ = j;
    }
    Matrix(int i, int j, double **matrix) {sizeI=i; sizeJ=j; pMatrix=matrix;}
    ~Matrix() {
        for (int i=0; i<GetSizeJ(); i++) delete [] pMatrix[i];
            delete [] pMatrix;
    }

    double Get(int i, int j) { return pMatrix[j][i]; }
    void Set(int i, int j, double value) { pMatrix[j][i]=value; }
    double* GetColumn(int j) { return pMatrix[j]; }

    int GetSizeI() {return sizeI; }
    int GetSizeJ() {return sizeJ; }

  Matrix* calculateMedian() {
    Matrix *res = new Matrix(GetSizeI(),1);
    double temp;
    // Check if there are columns to process.
    if (GetSizeJ()!=0) {
      for (int i=0;i<GetSizeI();i++) {
        temp = 0;

        for (int j=0;j<GetSizeJ();j++) {
          temp += Get(i,j);
        }
        temp = temp / GetSizeJ();

        res->Set(i,0,temp);
      }
    } else {
      for (int i=0;i<GetSizeI();i++) res->Set(i,0,0);
    }
    return (res);
  }

    friend std::ostream& operator<<(std::ostream& s, Matrix& p) {
        s << "Matrix(" << p.GetSizeI() << "," << p.GetSizeJ() << ") \n";
        for (int i=0;i<p.GetSizeI();i++) {
            for (int j=0;j<p.GetSizeJ();j++) {
                s << "|" << p.Get(i,j);
            }
            s << "|\n";
        }
        return s;
    }

  /*
  friend Matrix operator-(Matrix& s, Matrix& p) {

    if ((s.GetSizeI()!=p.GetSizeI())||(s.GetSizeJ()!=p.GetSizeJ())) throw new exception;

    Matrix r(p.GetSizeI(),p.GetSizeJ());
    double value;

    for (int i=0;i<s.GetSizeI();i++) {
      for (int j=0;j<s.GetSizeJ();j++) {
        value = s.Get(i,j) - p.Get(i,j);
        r.Set(i,j,value);
      }
    }

    return r;
  }
  */
};

/*!
 * \class Model
 * \brief Abstract class of Model
 */
class Model {
public:
private:
};

/*!
 * \class DataSource
 * \brief Abstract class of DataSource
 */
class DataSource {
public:
    //virtual double** ExtractFeatures (RawData &rw) =0;
    virtual Matrix* ExtractSampleFeatures () =0;
    virtual double* ExtractFeatures (Ptr<Packet> packet) =0;  //FIXME: Useless (was there before)
    virtual void UpdateDataBase(double *features) =0;
    virtual uint16_t GetNumberFeatures() {return m_num_features;};
    virtual ~DataSource();

protected:
    uint16_t m_num_features;
};

/*!
 * \class Detector
* \brief Abstract class of Detector
 */
class Detector {
public:
    //1st constructor
    Detector(DataSource *ds, Model *m, uint64_t p_init, uint64_t p_train_interval, bool p_only_one_train):
    m_datasource(ds), m_model(m), m_init_period(p_init), m_train_interval(p_train_interval), m_only_one_train(p_only_one_train) {};

    //2nd constructor
    Detector(DataSource *ds, uint64_t p_init, uint64_t p_train_interval, bool p_only_one_train):
    m_datasource(ds), m_model(0), m_init_period(p_init), m_train_interval(p_train_interval), m_only_one_train(p_only_one_train) {};

    virtual ~Detector();
    virtual Result Evaluate () =0;
    virtual void Train (bool anomaly) =0; // Executed each time new data is available.
    //virtual void IncrementalTrain (Ptr<Packet> packet, bool attack) =0;
    void registerFeatureExtractor(DataSource *datasource) {m_datasource = datasource;};

    //NOTE Disable since we doesn't send packets
    //virtual void SetPacket (Ptr<Packet> packet) =0;
    //virtual void UnSetPacket()=0;


protected:
    DataSource *m_datasource;
    Model *m_model;
    uint64_t m_init_period;
    uint64_t m_train_interval;
    bool m_only_one_train;

};

/*!
 * \class EuclideanModel
 * \brief Euclidian model for the 1st algorithm for AD
 * Inherit from Model 
 */
class EuclideanModel: public Model {

public:
  EuclideanModel() : max(0), min(0), threshold(0), avg_dst(0), dvt(0), average(0), num_elements(0), num_avg_elements(0), num_thr_elements(0) {
          NS_LOG_INFO("Creating model: " << this);
  }

  ~EuclideanModel() {
    NS_LOG_INFO("Removing model: " << this);
    if (average!=0) delete average;
    if (min!=0) delete [] min;
    if (max!=0) delete [] max;
    average = 0;
    min = 0;
    max = 0;
  }

  void SetMax(double *p_max) {
    max = p_max;
  }

  double* GetMax() {
    return max;
  }

  void SetMin(double *p_min) {
    min = p_min;
  }

  double* GetMin() {
    return min;
  }

  void SetThreshold(double p_thres) {
    threshold = p_thres;
  }

  double GetThreshold() {
    return threshold;
  }

  void SetAvgDst(double p_avg_dst) {
    avg_dst = p_avg_dst;
  }

  double GetAvgDst() {
    return avg_dst;
  }

  void SetDvt(double p_dvt) {
    dvt = p_dvt;
  }

  double GetDvt() {
    return dvt;
  }

  void SetAverage(Matrix *p_avg) {
    NS_LOG_INFO("Model " << this << " assigning matrix " << p_avg);
    average = p_avg;
  }

  Matrix* GetAverage() {
    return average;
  }

  void SetNumElements(int p_num_elements) {
    num_elements = p_num_elements;
  }

  int GetNumElements() {
    return num_elements;
  }

  void IncNumAvgElements() {
    num_avg_elements++;
  }

  int GetNumAvgElements() {
    return num_avg_elements;
  }

  void IncNumThrElements() {
    num_thr_elements++;
  }

  int GetNumThrElements() {
    return num_thr_elements;
  }

  void Reset() {
    num_elements = 0;
    num_avg_elements = 0;
    num_thr_elements = 0;
    threshold = 0;
    avg_dst = 0;
    dvt = 0;

    double* avg = average->GetColumn(1);

    for (int i = 0; i < average->GetSizeI(); i++) {
      avg[i] = 0;
      max[i] = 0;
      min[i] = 0;
    }
  }
    
  friend std::ostream& operator<<(std::ostream& s, EuclideanModel& p) {

    double* min = p.GetMin();
    double* max = p.GetMax();
    double* avg = p.GetAverage()->GetColumn(0);

    s << "Normality Model(T.Elem: " << p.GetNumElements() << ", Avg.Elem: " << p.GetNumAvgElements()
        << ", Thr.Elem: " << p.GetNumThrElements() << ", Threshold: " << p.GetThreshold() << ") \n\r";
    s << "Min: ";
    if (min!=0) {
      for (int i=0; i < p.GetAverage()->GetSizeI(); i++) {
        s << min[i] << ",";
      }
    }
    s << "\n\rMax: ";
    if (max!=0) {
      for (int i=0; i < p.GetAverage()->GetSizeI(); i++) {
       s << max[i] << ", ";
      }
    }
    s << "\n\rAverage: ";
        for (int i=0; i < p.GetAverage()->GetSizeI(); i++) {
         s << avg[i] << ", ";
        }
    s << "\n\r";
     return s;
   }

private:
  double *max, *min;
  double threshold;
  double avg_dst, dvt;
  Matrix *average;
  int num_elements;
  int num_avg_elements;
  int num_thr_elements;
    
};

/*!
 * \class DetectorEuclidean
 * \brief Euclidian detector used for the 1st algorithm for AD
 * Inherit from Detector
 */
class DetectorEuclidean : public Detector {
public:
    //TODO Understand why 2 constructors (x/o Model *m or not)
    //1st constructor w Model *m
    DetectorEuclidean(DataSource *ds, Ptr<Node> p_node, Model *m, uint64_t p_init, uint64_t p_train_interval, bool p_only_one_train, bool p_norm_data, bool p_time_ref):
    Detector(ds,m,p_init,p_train_interval,p_only_one_train), m_node(p_node) {};

    //2nd constructor w/o Model *m
    DetectorEuclidean(DataSource *ds, Ptr<Node> p_node, uint64_t p_init, uint64_t p_train_interval, bool p_only_one_train, bool p_norm_data, bool p_time_ref);

    ~DetectorEuclidean();

    Result Evaluate ();
    void Train (bool anomaly);
    //NOTE Disable since we doesn't send packets (see class Detector)
    //void SetPacket(Ptr<Packet> packet);
    //void UnSetPacket();
    //

    //void IncrementalTrain (Ptr<Packet> packet, bool attack);
    //double GetVectorNorm (const double *in_vector, int size);
    //void NormaliseVector (double *in_vector, int size);

    // Logging methods.
    void OpenLogFiles();
    void CloseLogFiles();
    void LogFeatures(double *features, int num_features);
    void LogNormFeatures(double *features, int num_features);
    void LogEvaluationResult(string packet_id, double distance, double threshold, double result, double* log_dst, int num_features, AttackType attack);
    void LogVectorProjections(double* projections, int size);

    // Training methods.
    void CreateNormalityModel();
    void UpdateThreshold(EuclideanModel *model, double* p_features, uint16_t p_num_feat, bool p_norm_data);
    void CalculateThreshold(EuclideanModel *model);
    double Distance(double *v_avg_feat, double *v_feat, int size, double* log_dst);
    void UpdateAverageMaxMin (EuclideanModel* p_inc_model, double* p_features, uint16_t p_num_feat);
    void CalculateAverage(EuclideanModel* p_inc_model, uint16_t p_num_feat, bool p_norm_data);
    void Normalise(double* max, double* min, double* vector, int num_elements);
    
private:
    // Variables
    Ptr<Node> m_node;
    ofstream f_log_norm_features;
    ofstream f_log_evaluation;
    ofstream f_log_projections;
    uint64_t num_packets, num_positives, num_not_eval;

    //Attack detection:
    double *m_packet_features;
    bool one_train;
    bool first_step_completed;
    bool m_norm_data;
    bool m_time_ref;
    EuclideanModel* m_inc_model;
    uint16_t m_num_feat;
};

// ---------- End of Classes -------------------------------------------------


// ---------- Main -----------------------------------------------------------
NS_LOG_COMPONENT_DEFINE ("ADSSimulator");

int 
main (int argc, char *argv[])
{
    NS_LOG_UNCOND ("ADS Simulator");

    // ---------- Simulation Variables -----------------------------------------
    // Change the variables and file names only in this block!
    //TODO Put variable that are needed for the simulation
    // ---------- End of Simulation Variables ----------------------------------


    // ---------- Network Setup ------------------------------------------------
    //TODO Setup the network (wave)
    // ---------- End of Network Setup -----------------------------------------

    // ---------- Allocate Node Positions --------------------------------------
    //TODO Setup the initialisation of the network
    // ---------- End of Allocate Node Positions -------------------------------

    // ---------- Use Mobibility Model for Node Positions ----------------------
    //TODO Affect Mobitlity to the node according to t
    // ---------- End of Mobibility Model --------------------------------------


    // ---------- Simulation Monitoring ----------------------------------------
    NS_LOG_INFO ("Configure Tracing.");
    //TODO: AsciiTraceHelper ascii

    Simulator::Run ();

    // Collect data there

    Simulator::Destroy ();
    // ---------- End of Simulation Monitoring ---------------------------------
    return 0;
}

// ---------- Function Definitions -------------------------------------------
// Detector
Detector::~Detector() {
    delete m_datasource;
    if (m_model!=0) delete m_model;
}

// DataSource
DataSource::~DataSource() {}

// DetectorEuclidean
/*!
 * Anomaly detector: Deviation from sum of squared differences
 * It creates the normality model in two steps. First calculates the averages
 * of the features and secondly calculates the threshold used for the distances.
 */
DetectorEuclidean::DetectorEuclidean(DataSource *ds, Ptr<Node> p_node, uint64_t p_init, uint64_t p_train_interval, bool p_only_one_train, bool p_norm_data, bool p_time_ref) : Detector(ds, p_init, p_train_interval, p_only_one_train) {
    //TODO: Review DetectorEuclidean::DetectorEuclidean constructor
    
    m_model = new EuclideanModel();
    m_node = p_node;
    one_train = false;
    first_step_completed = false;
    num_packets = 0;
    m_norm_data = p_norm_data;
    m_time_ref = p_time_ref;
    num_not_eval = 0;
    m_num_feat = ds->GetNumberFeatures();
    m_train_interval = p_train_interval;
    NS_LOG_DEBUG("Created model: " << m_model);
    //TODO: Create a loggging system for gather results and data
    //if (DEBUG_MODE) OpenLogFiles();

    // Initialise incremental model
    double* min = new double[m_num_feat];
    double* max = new double[m_num_feat];
    Matrix* avg = new Matrix(m_num_feat,1);
    m_inc_model = new EuclideanModel();
    m_inc_model->SetMax(max);
    m_inc_model->SetMin(min);
    m_inc_model->SetAverage(avg);
}

/*!
 * Destructor for DetectorEuclidean class
 */
DetectorEuclidean::~DetectorEuclidean() {
    //if (DEBUG_MODE) CloseLogFiles();
    if (m_inc_model!=0) delete m_inc_model;
  }

/*!
 * Function that evaluates the current observation against the normality model.
 * 
 * Returns: 0 - Normal evaluation
 *          1 - Abnormal evaluation
 */
Result DetectorEuclidean::Evaluate() {
    //TODO Do Result DetectorEuclidean::Evaluate (see ads-chisqr.cc)
    //TODO: Review DetectorEuclidean::Evaluate
    Result res;
    res.decision = false;
    res.believe = 2;
    res.attack = (AttackType) AT_MAXNUM;    
    return res;
}

/*!
 * Method that coordinates the training of the normality model.
 * 
 * Behaviour: Creates only one model after an initial period collecting data.
 */
void DetectorEuclidean::Train(bool anomaly) {
//TODO: Review DetectorEuclidean::Train
      
    num_packets++;

    if (!anomaly) {

      // Select source depending on the parameter units selected: time or packet number.
      uint64_t reference;
      if (m_time_ref) {

        reference = Simulator::Now().GetTimeStep();

        //if (reference > m_init_period) {

        if (((reference - m_init_period) < m_train_interval / 2)) {
          // Update average, max and min.
          UpdateAverageMaxMin(m_inc_model, m_packet_features, m_num_feat);

        }

        if (((reference - m_init_period) >= (m_train_interval / 2))&&((reference - m_init_period) < m_train_interval)) {

          if (!first_step_completed) {
            CalculateAverage(m_inc_model, m_num_feat, m_norm_data);

            first_step_completed = true;
          }

          UpdateThreshold(m_inc_model, m_packet_features, m_num_feat, m_norm_data);

        }
        // Check if the init period has been ended.
        // NS_LOG_DEBUG("Time ref: reference: " << reference << " init: " << m_init_period << " one train: " << m_only_one_train << " node " << m_node->GetId());
        if (reference >= (m_init_period + m_train_interval)) {

          //if ((m_only_one_train&&!one_train)||(!m_only_one_train)) {
          if (!one_train) {

           // Calculate threshold.
           CalculateThreshold(m_inc_model);

           // Apply new model.
           m_model = m_inc_model;
           m_inc_model = 0;
           one_train = true;

          }

        } else {
          num_not_eval++;
        }

        //}
      } else {
/*
        reference = num_packets;
        // Check if the init period has been ended.
        // NS_LOG_DEBUG("numpkts: " << num_packets << " init: " << m_init_period << " one train: " << m_only_one_train << " node " << m_node->GetId());
        if (reference >= m_init_period) {

          if ((m_only_one_train&&!one_train)||(!m_only_one_train)) {

            // Create the normality model after the init period.
            if (((reference - m_init_period)%m_train_interval) == 0) {
              CreateNormalityModel();
              one_train = true;
            }
          }

        } else {
          num_not_eval++;
        }
*/
      }

    }
  }

/*!
 * Calculate the threshold used to evaluate the distance.
 */
void DetectorEuclidean::CalculateThreshold(EuclideanModel *model) {
   //TODO: Review DetectorEuclidean::CalculateThreshold
 
    double avg_dst = model->GetAvgDst();
    double dvt = model->GetDvt();
    int num_elements = model->GetNumThrElements();

    // Calculate final deviation.
    dvt = sqrt(dvt / num_elements);
    model->SetDvt(dvt);

    // Calculate threshold.
    model->SetThreshold(avg_dst + dvt*3);

    // Update total number of elements used for training.
    model->SetNumElements(model->GetNumAvgElements() + num_elements);

    // Print debug info regarding the model.
    NS_LOG_DEBUG("Normality model Node " << m_node->GetId() << " Packet " << num_packets);
    NS_LOG_DEBUG((*model));
}

/*!
 * Calculate the threshold used to evaluate the distance.
 */
void DetectorEuclidean::UpdateThreshold(EuclideanModel *model, double* p_features, uint16_t p_num_feat, bool p_norm_data) {
//TODO: Review DetectorEuclidean::UpdateThreshold

    /*
     * To-Do: The algorithm calculates the variance, hence it should be stored
     *        in a variable called Variance and not Dvt as it is now.
     */

    double *average = model->GetAverage()->GetColumn(0);
    double *vector = p_features;
    double dst;
    double avg_dst = model->GetAvgDst();
    double dvt = model->GetDvt();
    double delta;

    // Normalise vector if needed
    if (p_norm_data) {
      vector = new double[p_num_feat];
      for (uint16_t i = 0; i < p_num_feat; i++) vector[i] = p_features[i];
      Normalise(model->GetMax(), model->GetMin(), vector, p_num_feat);
    }

    // Calculate Chi-Sqr distance of the current observation.
    dst = Distance(average, vector, p_num_feat, 0);

    // Update average and deviation (Knuth algorithm, see
    // Wikipedia http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance)
    model->IncNumThrElements();
    delta = dst - avg_dst;
    avg_dst = avg_dst + delta / model->GetNumThrElements();
    dvt = dvt + delta * (dst - avg_dst);

    // Update model values.
    model->SetAvgDst(avg_dst);
    model->SetDvt(dvt);

    // Delete created vector if normalisation.
    if (p_norm_data) delete[] vector;

  }

/*!
 * Calculate distance between two vectors based on sum of squared differences.
 */
double DetectorEuclidean::Distance(double *v_1, double *v_2, int size, double *log_dst) {
//TODO: Review DetectorEuclidean::Distance

    double result = 0;
    double dst = 0;

    for (int i=0; i<size; i++) {
    //NS_LOG_DEBUG("ChiSquare: Feature:  " << features[i] << "  Average: " << average[i] << "  subtract: " << (features[i] - average[i]) << " part: " << (pow((features[i] - average[i]),2)/average[i]) << " presult: " << result);
    // NS_LOG_WARN("ChiSquare: Feature[" << i << "]: " << features[i] << "  Average: " << average[i]);
    // dst = pow((v_feat[i] - v_avg_feat[i]),2)/v_avg_feat[i];  <-- Before making average absolute in the denominator
    dst = v_1[i]-v_2[i];
    //if (log_dst != 0) log_dst[i] = dst; // Original log
    if (log_dst != NULL) log_dst[i] = dst;  // Log for Tiziano
    result += dst*dst; //multiplication instead of pow()
    }

    // TESTING REAL EUCLIDEAN!!!!! (It works worse. Long distances are more relevant without SQRT)
    result = sqrt(result);

    return result;

   }

/*!
 * Method to normalise a vector.
 */
void DetectorEuclidean::Normalise(double* max, double* min, double* vector, int num_elements) {
//TODO: Review DetectorEuclidean::Normalise
    for (int i = 0; i < num_elements; i ++) {
        if ((max[i] - min[i]) == 0) {
            vector[i] = 0;
        } else {
            //if (min[i]<=0) vector[i] = (vector[i] + min[i])/(max[i] - min[i]);
            //if (min[i]>0) vector[i] = (vector[i] - min[i])/(max[i] - min[i]);
            vector[i] = (vector[i] - min[i])/(max[i] - min[i]);
            // if (features[i] < 0) features[i] = 0;
            // if (features[i] > 1) features[i] = 1;
        }
    }

}

/*!
 * Method to update the temporal model with the new observation.
 *
 * Notice that the average vector is only incremented instead of averaging it.
 * This is because of performance reasons and for getting maximum precision.
 */
void DetectorEuclidean::UpdateAverageMaxMin (EuclideanModel* p_inc_model, double* p_features, uint16_t p_num_feat) {
//TODO: Review DetectorEuclidean::UpdateAverageMaxMin
    double* avg = p_inc_model->GetAverage()->GetColumn(0);
    double* max = p_inc_model->GetMax();
    double* min = p_inc_model->GetMin();

    // Add the current observation and check max and min values.
    if (p_inc_model->GetNumAvgElements() == 0) {
        for (uint16_t i = 0; i < p_num_feat; i++) {
            avg[i] = p_features[i];
            max[i] = p_features[i];
            min[i] = p_features[i];
        }
    } else {
        for (uint16_t i = 0; i < p_num_feat; i++) {
            avg[i] = avg[i] + p_features[i];
            if (p_features[i] > max[i]) max[i] = p_features[i];
            if (p_features[i] < min[i]) min[i] = p_features[i];
        }
    }

    // Increase the number of observations registered.
    p_inc_model->IncNumAvgElements();


}

/*!
 * Method that calculates final average vector.
 */
void DetectorEuclidean::CalculateAverage(EuclideanModel* p_inc_model, uint16_t p_num_feat, bool p_norm_data) {
    //TODO: Review DetectorEuclidean::CalculateAverage
    double* avg = p_inc_model->GetAverage()->GetColumn(0);
    double* max = p_inc_model->GetMax();
    double* min = p_inc_model->GetMin();

    int num_elements = p_inc_model->GetNumAvgElements();

    for (uint16_t i = 0; i < p_num_feat; i ++) {
        avg[i] = avg[i] / num_elements;
    }

    if (p_norm_data) {
        Normalise(max, min, avg, p_num_feat);
    }
}

// ---------- End of Function Definitions ------------------------------------









