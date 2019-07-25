/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pprfinder;

/**
 *
 * @author ian
 */
public class RunningStat
    {
        private int m_n;
        private double m_oldM, m_newM, m_oldS, m_newS;

        void Clear()
        {
            m_n = 0;
        }

        void Push(double x)
        {
            m_n++;

            // See Knuth TAOCP vol 2, 3rd edition, page 232
            if (m_n == 1)
            {
                m_oldM = m_newM = x;
                m_oldS = 0.0;
            }
            else
            {
                m_newM = m_oldM + (x - m_oldM)/m_n;
                m_newS = m_oldS + (x - m_oldM)*(x - m_newM);
    
                // set up for next iteration
                m_oldM = m_newM; 
                m_oldS = m_newS;
            }
        }

        public int NumDataValues()
        {
            return m_n;
        }

        public double Mean()
        {
            return (m_n > 0) ? m_newM : 0.0;
        }

        public double Variance()
        {
            return ( (m_n > 1) ? m_newS/(m_n - 1) : 0.0 );
        }

        public double StandardDeviation()
        {
            return Math.sqrt( Variance() );
        }
        
    };
