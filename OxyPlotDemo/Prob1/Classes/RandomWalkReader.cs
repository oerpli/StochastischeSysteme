using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using OxyPlot;
using StochSyst.ViewModels;
using System.Globalization;

namespace StochSyst {
    public abstract class RandomWalkReader {
        protected CultureInfo culture = new CultureInfo("en");

        protected int omit;
        protected int max;
        protected string filename;
        protected MainWindowModel window;
        protected RandomWalkPlot plot;
        protected RandomWalkReader(MainWindowModel mw, string filename, string description, int omit, int max) {
            this.max = max;
            this.filename = filename;
            this.window = mw;
            this.omit = omit;
            plot = new RandomWalkPlot(window.Plot1, "random walk from file: " + filename, description + (omit > 1 ? " (every " + omit + "th element plotted)" : "."));
            ReadFile();
        }

        private void ReadFile() {
            System.IO.StreamReader file = new System.IO.StreamReader(filename);
            this.PlotFile(file);
            file.Close();
        }

        private void PlotFile(System.IO.StreamReader file) {
            string line;
            int count = 0;
            while((line = file.ReadLine()) != null && count <= max) {
                if(count == 0) {
                    InitPlot("", OxyColors.Red);
                }
                count++;
                if(count % omit == 0) {
                    PlotLine(line, count);
                }
            }
            //AddDistPlot();
            //DrawHistogram();
        }


        public abstract void PlotLine(string line, int lineIndex);
        public abstract void InitPlot(string firstline, OxyColor? color);
        public abstract void AddDistPlot();


        private int[,] hist = new int[40, 1500];


        protected void HistogramData(int time, double distance) {
            int dist = (int)distance;// ((int)(distance * 6) / 3) / 2;
            hist[time / 10, dist]++;
        }
        private void DrawHistogram() {
            //plot.AddLine("Hist10", OxyColors.Black);
            //for(int i = 0; i < hist10.Length; i++) {
            //    if(hist10[i]>0)plot.AddPoint(i,hist10[i]);
            //}
            for(int j = 1; j < 35; j+=5) {
                byte c =(byte)( j * 5);
                plot.AddLine("Hist (N=" + j*10 + ")", OxyColor.FromRgb(c,c,c));
                for(int i = 0; i < 1500; i++) {
                    if(hist[j, i] > 0) plot.AddPoint(i, hist[j, i]);
                }
            }
            int x = 34;
            plot.AddLine("Hist" + x * 10, OxyColors.Red);
            for(int i = 0; i < 1500; i++) {
                if(hist[x, i] > 0) plot.AddPoint(i, hist[x, i]);
            }
        }

    }


    //public class RandomWalkXT : RandomWalkReader
    //{
    //    private RandomWalkPlot[] plots;
    //    public RandomWalkXT(MainWindowModel mw, string filename, string description, int omit) : base(mw, filename, description, omit) { }

    //    public override void InitPlot(string filename, string firstline)
    //    {
    //        String[] values = firstline.Split(' ');
    //        plots = new RandomWalkPlot[values.Length];
    //        for (int i = 0; i < values.Length; i++)
    //        {
    //            plots[i] = rwplot.AddLine(filename + "_" + i, OxyColors.MediumBlue);
    //        }

    //    }
    //    //public override void PlotLine(string line, int lineIndex)
    //    //{
    //    //    String[] values = line.Split(' ');
    //    //    double[] doubles = new double[values.Length];
    //    //    for (int i = 0; i < values.Length; i++)
    //    //    {
    //    //        Console.WriteLine(line);
    //    //        doubles[i] = double.Parse(values[i], culture);
    //    //        plots[i].AddPoint(lineIndex, doubles[i]);
    //    //    }
    //    //}
    //}

    public class RandomWalk2D : RandomWalkReader {
        private Random rng = new Random();
        //private RandomWalkPlot plot;
        private bool omitnext = false;
        private int time = 0;
        private int plots = 0;
        private List<Double> avgdist = new List<Double>();

        public RandomWalk2D(MainWindowModel mw, string filename, string description, int omit, int max)
            : base(mw, filename, description, omit, max) { }

        public override void InitPlot(string Name, OxyColor? color) {
            plot = plot.AddLine(Name, color ?? OxyColors.DodgerBlue);
            plots++;
        }
        public override void PlotLine(string line, int lineIndex) {
            String[] values = line.Split(' ');
            if(omitnext || values.Length != 2) {
                InitPlot("", OxyColor.FromUInt32((uint)rng.Next()));
                time = 0;
                omitnext = !omitnext;
                return;
            }
            double[] doubles = new double[values.Length];


            for(int i = 0; i < values.Length; i++) {
                doubles[i] = double.Parse(values[i], culture);
            }


            double dist = Math.Sqrt(Math.Pow(doubles[0], 2) + Math.Pow(doubles[1], 2));

            //Split single walk 
            //if(dist < 0.05 && time > 1000) {
            //    InitPlot("", OxyColor.FromUInt32((uint)rng.Next()));
            //    time = 0;
            //}

            if(avgdist.Count <= time)
                avgdist.Add(0);
            avgdist[time] += dist;

            plot.AddPoint(time, dist);
            //plot.AddPoint(doubles[0], doubles[1]);
            if(time % 10 == 0) {
                HistogramData(time, dist);
            }
            time += omit;
        }
        public override void AddDistPlot() {
            InitPlot("E(X²+Y²)", OxyColors.Red);
            int listlength = avgdist.Count;
            for(int i = 0; i < listlength; i += omit) {
                plot.AddPoint(i, avgdist[i] / plots);
            }
        }
    }
}

