/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
#ifndef AD_DETECTOR_SIMULATOR_H
#define AD_DETECTOR_SIMULATOR_H

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
using namespace std;
namespace ns3 {
    // ---------- Prototypes -----------------------------------------------------
    class Matrix;
    class Model;
    class DataSource;
    class Detector;
    class EuclideanModel;
    class DetectorEuclidean;
    double vDotProduct (const double *a_vector,const double *b_vector, int size);
    double GetVectorNorm (const double *in_vector,int size);
    void NormaliseVector (double *in_vector,int size);
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
          //NS_LOG_INFO("Creating model: " << this);
      }

      ~EuclideanModel() {
        //NS_LOG_INFO("Removing model: " << this);
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
        //NS_LOG_INFO("Model " << this << " assigning matrix " << p_avg);
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
        double avg_dst;     //!< Average Distance
        double dvt;         //!< Deviation 
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
        // 2 constructors available according to the situation:
        //either we have a Model *m already existing for to the new detector
        //either we create the model in the constructor if not given as parameter
        //1st constructor w Model *m
        DetectorEuclidean(DataSource *ds, Ptr<Node> p_node, Model *m, uint64_t p_init, uint64_t p_train_interval, bool p_only_one_train, bool p_norm_data, bool p_time_ref):
        Detector(ds,m,p_init,p_train_interval,p_only_one_train), m_node(p_node) {};

        //2nd constructor w/o Model
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
    
    
} // namespace ns3

#endif /* AD_DETECTOR_SIMULATOR_H */

