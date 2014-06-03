using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using OxyPlot.Series;
using OxyPlot.Axes;
using OxyPlot;
using StochSyst.ViewModels;

namespace StochSyst {

    class Sampler {
        static protected Random rng = new Random();
        private const int N = 2000;
        private double[,] samples = new double[2, N];
        private MainWindowModel window;
        private RandomWalkPlot rwplot;

        public Sampler(MainWindowModel mw) {
            this.window = mw;
            rwplot = new RandomWalkPlot(window.Plot1, "Sampling Gauss/Cauchy ", "");
            for(int i = 0; i < N; i++) {
                samples[0, i] = cauchyRandom();
                samples[1, i] = gaussRandom();
            }
            Console.WriteLine(samples[0, 2]);
            rwplot.AddLine("Cauchy", OxyColors.Blue);
            for(int i = 0; i < N; i++) {
                double mean= 0;
                for(int j = 0; j < i; j++) { 
                    mean = mean + samples[0,j];
                }
                rwplot.AddPoint(i,mean/i);
            }
            rwplot.AddLine("Gauss", OxyColors.Red);
            for(int i = 0; i < N; i++) {
                double mean = 0;
                for(int j = 0; j < i; j++) {
                    mean = mean + samples[1, j];
                }
                rwplot.AddPoint(i, mean / i);
            }
        }

        private double cauchyRandom() {
            return (Math.Tan((rng.NextDouble() - 0.5) * Math.PI));
        }
        private double gaussRandom() {
            double u1 = rng.NextDouble();
            double u2 = rng.NextDouble();
            return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2.0 * Math.PI * u2);
        }
    }
}
