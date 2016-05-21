#ifndef C_THREAD_DELTAU_MOVE_VERTEX_H
#define C_THREAD_DELTAU_MOVE_VERTEX_H

#include <C_thread.h>
#include <C_energy.h>

class C_energy;

class C_thread_deltaU_move_vertex : public C_thread
{
    //friend class C_energy;

    public:
        /** Default constructor */
        C_thread_deltaU_move_vertex(C_energy* U);
        /** Default destructor */
        virtual ~C_thread_deltaU_move_vertex();

        virtual void execute();

        void resetAttributes();

        Point* u; //!< Member variable "u"
        Point* v; //!< Member variable "v"
        Point* v_new; //!< Member variable "v_new" this point is not part of the graph, it does not have all energy caluculated
        vector<long> idxPair; //!< Member variable "idxPair"
        vector<double> dataPair; //!< Member variable "dataPair"
        double dU_data; //!< Member variable "dU_data"
        //double dU_topo; //!< Member variable "dU_topo"
        unsigned long m_kind;

    protected:
        C_energy* m_U;

};

#endif // C_THREAD_DELTAU_MOVE_VERTEX_H










//        /** Access u
//         * \return The current value of u
//         */
//        Point* Getu() { return u; }
//        /** Set u
//         * \param val New value to set
//         */
//        void Setu(Point* val) { u = val; }
//        /** Access v
//         * \return The current value of v
//         */
//        Point* Getv() { return v; }
//        /** Set v
//         * \param val New value to set
//         */
//        void Setv(Point* val) { v = val; }
//        /** Access v_new
//         * \return The current value of v_new
//         */
//        double* Getv_new() { return v_new; }
//        /** Set v_new
//         * \param val New value to set
//         */
//        void Setv_new(double* val) { v_new = val; }
//        /** Access idxPair
//         * \return The current value of idxPair
//         */
//        vector<unsigned long> GetidxPair() { return idxPair; }
//        /** Set idxPair
//         * \param val New value to set
//         */
//        void SetidxPair(vector<unsigned long> val) { idxPair = val; }
//        /** Access dataPair
//         * \return The current value of dataPair
//         */
//        vector<double> GetdataPair() { return dataPair; }
//        /** Set dataPair
//         * \param val New value to set
//         */
//        void SetdataPair(vector<double> val) { dataPair = val; }
//        /** Access dU_data
//         * \return The current value of dU_data
//         */
//        double GetdU_data() { return dU_data; }
//        /** Set dU_data
//         * \param val New value to set
//         */
//        void SetdU_data(double val) { dU_data = val; }
//        /** Access dU_topo
//         * \return The current value of dU_topo
//         */
//        double GetdU_topo() { return dU_topo; }
//        /** Set dU_topo
//         * \param val New value to set
//         */
//        void SetdU_topo(double val) { dU_topo = val; }
