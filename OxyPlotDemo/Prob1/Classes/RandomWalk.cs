using System.Globalization;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using OxyPlot.Series;
using OxyPlot.Axes;
using OxyPlot;
using StochSyst.ViewModels;
using MathNet.Numerics.Signals;
using MathNet.Numerics.IntegralTransforms;
using System.Numerics;


namespace StochSyst {
    public class RandomWalk {
        static protected Random rng = new Random();

        private MainWindowModel window;
        private RandomWalkPlot rwplot;
        protected const int Traj = 2 * 1000 * 1000;
        protected const int Steps = 100;
        private double[] mean = new double[101];
        protected const int reflPos = 30;
        protected const int classesHist = 200;
        protected int[,] positions = new int[classesHist, 11];


        //private double[,] hist = new double[4, 800];
        public RandomWalk(MainWindowModel mw) {
            this.window = mw;
            //rwplot = new RandomWalkPlot(window.Plot1, "Diffusion", " with reflective boundaries, analytical solution");
            //+ (Traj / 1000) + "k trajectories ");
            //TakeWalk(new WalkerCauchy(), OxyColors.Red);
            //TakeWalk(new WalkerBinomial(), OxyColors.RoyalBlue);

            //rwplot.SetAxisTitles("Position x", "p(x,t)");
            //OxyColor[] cols = { OxyColors.Blue, OxyColors.Red, OxyColors.Green, OxyColors.Black };
            //String[] titles = { "Reflect", "Truncate", "Discard", "Redraw" };

            //rwplot.AddLine(titles[0], OxyColors.DodgerBlue);
            //byte handling = 0;
            //double reflpos = 4.9;


            //2.2a
            //for(double reflpos = 1; reflpos <= reflPos; ) {
            //Console.WriteLine(reflpos + " / " + reflPos);
            //byte collol = (byte)(reflpos * 8);
            //rwplot.AddLine("Corresponding solution", OxyColors.Black);
            //positions = new int[classesHist, 11];
            //TakeWalk(new WalkerGaussianReflective(handling, -10, 10), cols[handling]);
            //rwplot.AddPoint(0, 0);
            //for(int step = 1; step <= 10; step++) {
            //    byte collol = (byte)(step * 20);
            //    rwplot.AddLine("",OxyColors.Black);
            //    for(int xcord = 0; xcord < classesHist; xcord++) {
            //        rwplot.AddPoint(xcord * scale + minX,
            //            reflektierteverteilung(0.5, xcord * scale + minX, step * step, -10, 10));
            //    }
            //    rwplot.AddLine("N=" + step * step, OxyColor.FromRgb(collol, (byte)(255 - collol), (byte)(255 - collol)));
            //    for(int xcord = 0; xcord < classesHist; xcord++) {
            //        rwplot.AddPoint(xcord * scale + minX + 0.5*scale, (1 / scale) * ((double)positions[xcord, step]) / Traj);
            //    }
            //}
            //if(reflpos < 4.95 && reflpos > 3.9) reflpos += 0.1;
            //else reflpos++;
            //}


            //2.1a
            //double minnX = -5;
            //double scalee = 0.01;
            //double[] temps = { 1, 10, 20, 30, 40, 50, 100 };
            //for(int i = 0; i < temps.Length; i++) {
            //    double T = temps[i];
            //    byte collol = (byte)(T * 30); OxyColor col = OxyColor.FromRgb(0, (byte)(255 - collol), (byte)(255 - collol));
            //    rwplot.AddLine("t=" + T, col);
            //    for(double xcord = minnX; xcord < -minnX; xcord += scalee) {
            //        rwplot.AddPoint(xcord,
            //            reflektierteverteilung(0.1, xcord, T, minnX, -minnX));
            //    }
            //}

            ////2.1b
            //rwplot.SetAxisTitles("Dt/L²", "<x²>/L²");
            //double scalee = 0.01;
            //double[] temps = { 1, 10, 20, 30, 40, 50, 100 };
            //double D = 2;
            //double L = 20;
            //double F = D / L / L;
            //double bla = reflektiertmsd(D, 1000000, L);
            //rwplot.AddLine("<x²> t = ∞", OxyColors.DodgerBlue);
            //for(double T = 0; T * F < 0.12; T += scalee) {
            //    rwplot.AddPoint(F * T, bla);
            //}
            //bla = reflektiertmsd(D, 1, L);
            //rwplot.AddLine("<x²> freies Teilchen", OxyColors.Red);
            //for(double T = 0; T * F < 0.12; T += scalee) {
            //    rwplot.AddPoint(F * T, bla*T);
            //}
            //rwplot.AddLine("<x²>", OxyColors.Black);
            //for(double T = 0; T * F < 0.12; T += scalee) {
            //    rwplot.AddPoint(F * T, reflektiertmsd(D, T, L));
            //}


            //4.1a
            //int yCount = 500000;
            //rwplot = new RandomWalkPlot(window.Plot1, "Sum of uniformly distributed values", yCount / 100000 + "E5 Samples");
            //rwplot.SetAxisTitles("n", "<y^m>");
            //int[] moments = { 4, 6 };
            //LineSeries[] plots = new LineSeries[moments.Length];
            //for(int i = 0; i < moments.Length; i++) {
            //    byte collol = (byte)(i * 40); OxyColor col = OxyColor.FromRgb(0, (byte)(255 - collol), (byte)(255 - collol));
            //    plots[i] = rwplot.AddLineM("Moment" + moments[i], col);
            //}
            //for(int n = 1; n <= 30; n++) {
            //    double mean = 0;
            //    double sumMoment = 0;

            //    //calc y-s
            //    double[] ys = new double[yCount];
            //    for(int i = 0; i < yCount; i++) {
            //        ys[i] = SumUniform(n);
            //        mean += ys[i];
            //    }
            //    mean /= yCount;

            //    //Calculate Moments
            //    for(int i = 0; i < moments.Length; i++) {
            //        for(int j = 0; j < yCount; j++) {
            //            sumMoment += Math.Pow(ys[j] - mean, moments[i]);
            //        }
            //        sumMoment /= yCount;
            //        rwplot.AddPoint(plots[i], n, sumMoment);
            //    }
            //}

            ////Asymptote
            //rwplot.AddLine("asympt. Moment 4", OxyColors.Black);
            //for(int n = 1; n <= 30; n++) {
            //    rwplot.AddPoint(n, 3);
            //} rwplot.AddLine("asympt. Moment 6", OxyColors.Black);

            //for(int n = 1; n <= 30; n++) {
            //    rwplot.AddPoint(n, 15);
            //}



            //4.1b
            //rwplot = new RandomWalkPlot(window.Plot1, "Sum of uniformly distributed values", "1E7 Samples");
            //rwplot.SetAxisTitles("y-6*n", "p(y+6*n)");
            //int[] ns = { 1, 2, 3, 5, 10 };
            ////int[] ns = {3};
            //var hist = new RectangleBarSeries[ns.Length];

            //double mean = 0;
            //for(int n = 0; n < ns.Length; n++) {
            //    int yCount = 2000000;
            //    double[] ys = new double[yCount];
            //    for(int i = 0; i < yCount; i++) {
            //        ys[i] = SumUniform(ns[n]);
            //        mean += ys[i];
            //    }
            //    mean = mean / yCount;


            //    hist[n] = new RectangleBarSeries();
            //    hist[n].Title = "Histogram for n = " + ns[n];
            //    int classN = 200;
            //    double intervall = 10;
            //    double isteps = 2 * intervall / classN;
            //    double[] classes = new double[classN];
            //    int counted = 0;
            //    for(int i = 0; i < yCount; i++) {
            //        for(int c = 0; c < classN; c++) {
            //            double classMin = -intervall + c * isteps;
            //            if(ys[i] >= classMin && ys[i] < classMin + isteps) {
            //                classes[c] += 1;
            //                counted++;
            //                break;
            //            }
            //        }
            //    }
            //    for(int c = 0; c < classN; c++) {
            //        double classMin = -intervall + c * isteps;
            //        if(classes[c] > 0)
            //            hist[n].Items.Add(new RectangleBarItem() { X0 = classMin + 6 * n, X1 = classMin + isteps + 6 * n, Y0 = 0, Y1 = classes[c] / yCount / isteps, });
            //    }

            //    window.Plot1.Series.Add(hist[n]);
            //    rwplot.AddLine("", OxyColors.Black);
            //    for(int c = 0; c < classN; c++) {
            //        double classMin = -intervall + c * isteps;
            //        if(classes[c] > 0)
            //            rwplot.AddPoint(classMin + isteps * 0.5 + 6 * n, normalverteilung(classMin + isteps * 0.5, 0, 1));
            //    }

            //}
            //rwplot.AddLine("Gaussian, µ=0, s = 1", OxyColors.Black);




            //4.2
            //    rwplot = new RandomWalkPlot(window.Plot1, "Dellago distributed values", "5E4 Samples");
            //    rwplot.SetAxisTitles("y-6*n", "p(y+6*n)");
            //    int[] ns = { 0 };
            //    //int[] ns = {3};
            //    var hist = new RectangleBarSeries[ns.Length];

            //    double mean = 0;
            //    for(int n = 0; n < ns.Length; n++) {
            //        int yCount = 500000;
            //        double[] ys = new double[yCount];
            //        for(int i = 0; i < yCount; i++) {
            //            ys[i] = DellagoRandom();
            //            mean += ys[i];
            //        }
            //        mean = mean / yCount;


            //        hist[n] = new RectangleBarSeries();
            //        hist[n].Title = "Histogram for n = " + ns[n];
            //        int classN = 200;
            //        double intervall = 4;
            //        double isteps = 2 * intervall / classN;
            //        double[] classes = new double[classN];
            //        int counted = 0;
            //        for(int i = 0; i < yCount; i++) {
            //            for(int c = 0; c < classN; c++) {
            //                double classMin = -intervall + c * isteps;
            //                if(ys[i] >= classMin && ys[i] < classMin + isteps) {
            //                    classes[c] += 1;
            //                    counted++;
            //                    break;
            //                }
            //            }
            //        }
            //        for(int c = 0; c < classN; c++) {
            //            double classMin = -intervall + c * isteps;
            //            if(classes[c] > 0)
            //                hist[n].Items.Add(new RectangleBarItem() { X0 = classMin, X1 = classMin + isteps, Y0 = 0, Y1 = classes[c] / yCount / isteps, });
            //        }
            //        window.Plot1.Series.Add(hist[n]);
            //    }
            //    rwplot.AddLine("Dellago", OxyColors.Black);
            //    double x = -5;
            //    while(x < 5) {
            //        rwplot.AddPoint(x, DellagoDistrib(x));
            //        x += 0.01;
            //    }
            //    x = -5;

            //    rwplot.AddLine(multi + "*Gaussian, µ=±1, s² = " + variance, OxyColors.Red);
            //    while(x < 5) {
            //        rwplot.AddPoint(x, DellagoMajor(x));
            //        //rwplot.AddPoint(x, normalverteilung(x, 0, 1));
            //        x += 0.01;
            //    }
            //    rwplot.AddLine(((double)rej / (double)acc).ToString() + "= Rej/Acc", OxyColors.White);

            //4.3
            //rwplot = new RandomWalkPlot(window.Plot1, "Stationary (?) process", "1E3 x 1E4 Samples");
            //CultureInfo culture = new CultureInfo("en");
            //string filename = "station.dat";
            //string line = "";
            //int histC = 0;
            //double timeM = 0;
            //int[] histT = { 2000, 2500, 3000, 3500 };
            //var hist = new RectangleBarSeries[histT.Length];
            //double max = 0;
            //double min = 0;
            ////<X(t)>
            //System.IO.StreamReader file = new System.IO.StreamReader(filename);

            //int xlol = 2;

            //bool xt = xlol == 0;
            //bool xhist = xlol == 1;
            //bool xauto = xlol == 2;
            //if(xt) rwplot.AddLine("<X(t)>", OxyColors.Black);
            //line = file.ReadLine();
            //int count = 1;
            //while(xt && (line = file.ReadLine()) != null) {
            //    rwplot.SetAxisTitles("t", "x(t)");
            //    string[] cols = line.Split(' ');
            //    timeM = 0;
            //    for(int col = 0; col < cols.Length; col++) {
            //        double x = double.Parse(cols[col], culture);
            //        timeM += x;
            //        if(x > max) max = x;
            //        if(x < min) min = x;
            //    }
            //    rwplot.AddPoint(count, timeM / cols.Length);
            //    count++;
            //}
            //file.Close();
            ////histogramme
            //count = 1;
            //file = new System.IO.StreamReader(filename);
            //line = file.ReadLine();
            //while(xhist && (line = file.ReadLine()) != null) {
            //    rwplot.SetAxisTitles("x", "p(x)");
            //    if(count == 2000 || count == 2500 || count == 3000 || count == 3500) {
            //        hist[histC] = new RectangleBarSeries();
            //        hist[histC].Title = "Histogram for n = " + count;
            //        int classN = 400;
            //        double intervall = 2;
            //        double isteps = 2 * intervall / classN;
            //        double[] classes = new double[classN];

            //        int classmean = 0;

            //        string[] cols = line.Split(' ');
            //        for(int col = 0; col < cols.Length; col++) {
            //            double x = double.Parse(cols[col], culture);
            //            for(int c = 0; c < classN; c++) {
            //                double classMin = -intervall + c * isteps;
            //                if(x >= classMin && x < (classMin + isteps)) {
            //                    classes[c] += 1;
            //                    classmean += c;
            //                    break;
            //                }
            //            }
            //        }
            //        Console.WriteLine(classmean);
            //        for(int c = 0; c < classN; c++) {
            //            double classMin = -intervall + c * isteps;
            //            if(classes[c] > 0)
            //                hist[histC].Items.Add(new RectangleBarItem() { X0 = classMin, X1 = classMin + isteps, Y0 = 0, Y1 = classes[c] / 1000 / isteps });
            //        }
            //        window.Plot1.Series.Add(hist[histC]);
            //        //histC++;
            //    }
            //    count++;
            //}
            //Console.WriteLine(count);
            //file.Close();

            ////Autocorrelation
            //OxyColor[] colors = { OxyColors.Blue, OxyColors.Red, OxyColors.Green, OxyColors.Black };
            //count = 0;
            //int autoC = 0;
            //file = new System.IO.StreamReader(filename);
            //if(xauto) rwplot.SetAxisTitles("a", "<X(t)X(t+a)");
            //while(xauto && (line = file.ReadLine()) != null) {
            //    const int corrlength = 100;
            //    double[,] vals = new double[corrlength, 1000];
            //    if(count == 2000 || count == 2500 || count == 3000 || count == 3500) {
            //        string[] cols = line.Split(' ');
            //        for(int col = 1; col < cols.Length; col++) {
            //            vals[0, col] = double.Parse(cols[col], culture);
            //        }
            //        for(int i = 1; i < corrlength; i++) {
            //            cols = file.ReadLine().Split(' ');
            //            for(int col = 0; col < cols.Length; col++) {
            //                vals[i, col] = double.Parse(cols[col], culture);
            //            }
            //        }//corrlength*1000 Values Stored



            //        byte collol = (byte)(70 * autoC++);
            //        rwplot.AddLine("<X(t)X(t+a)> - t=" + count, OxyColor.FromArgb(255, 0, (byte)(255 - collol), (byte)(255 - collol)));
            //        //calc mean & var
            //        double[] meanI = new double[1000];
            //        double[] varI = new double[1000];
            //        for(int col = 0; col < 1000; col++) {
            //            for(int i = 0; i < corrlength; i++) {
            //                meanI[col] += vals[i, col];
            //            }
            //        }
            //        for(int col = 0; col < 1000; col++) {
            //            meanI[col] /= corrlength;
            //        }
            //        for(int col = 0; col < 1000; col++) {
            //            for(int i = 0; i < corrlength; i++) {
            //                varI[col] += Math.Pow(vals[i, col] - meanI[col], 2);
            //            }
            //        }

            //        //Calculate Autocorrelation
            //        double[,] autocorr = new double[corrlength, 1000];
            //        for(int col = 0; col < 1000; col++) {
            //            for(int i = 0; i < corrlength; i++) {
            //                for(int j = 0; j < corrlength - i; j++) {
            //                    autocorr[j, col] += (vals[j + i, col] - meanI[col]) * (vals[i, col] - meanI[col]);
            //                    //autocorr[j, col] += (vals[j + i, col] * vals[i, col]);
            //                }
            //            }
            //        }
            //        for(int col = 0; col < 1000; col++) {
            //            for(int i = 0; i < corrlength; i++) {
            //                autocorr[i, col] /= varI[col];
            //            }
            //        }
            //        for(int col = 1; col < 1000; col++) {
            //            for(int i = 0; i < corrlength; i++) {
            //                autocorr[i, 0] += autocorr[i, col];
            //            }
            //        }
            //        for(int i = 0; i < corrlength; i++) {
            //            rwplot.AddPoint(i, autocorr[i, 0] / 1000);
            //        }
            //        //Finished
            //        count += corrlength - 1;
            //    }
            //    count++;
            //}
            //file.Close();



            //4.4
            //int corrlength = 1000;

            //rwplot = new RandomWalkPlot(window.Plot1, "Harmonic trap", "1E6 Samples, max correlation length = " + corrlength);
            //CultureInfo culture = new CultureInfo("en");
            //string filename = "harmonic.dat";
            //string line = "";
            //double timeM = 0;
            //int[] histT = { 2000, 2500, 3000, 3500 };
            //var hist = new RectangleBarSeries[histT.Length];
            //System.IO.StreamReader file = new System.IO.StreamReader(filename);


            //int xlol = 1;

            //bool xt = xlol == 0;
            //bool xauto = xlol == 1;

            ////<X(t)>
            //if(xt) rwplot.AddLine("<X(t)>", OxyColors.Black);
            //line = file.ReadLine();
            //int count = 1;
            //while(xt && (line = file.ReadLine()) != null) {
            //    rwplot.SetAxisTitles("t", "x(t)");
            //    double x = double.Parse(line, culture);
            //    timeM += x;
            //    rwplot.AddPoint(count, x);
            //    count++;
            //}
            //file.Close();


            ////Autocorrelation
            //OxyColor[] colors = { OxyColors.Blue, OxyColors.Red, OxyColors.Green, OxyColors.Black };
            //count = 0;
            //file = new System.IO.StreamReader(filename);
            //if(xauto) {
            //    rwplot.SetAxisTitles("a", "<X(t)X(t+a)");

            //    double[] vals = new double[100000];
            //    corrlength = 2000;
            //    for(int i = 1; i < vals.Length; i++) {
            //        vals[i] = double.Parse(file.ReadLine(), culture);
            //    }//all values stored.
            //    file.Close();

            //    //Calculate Autocorrelation O(n^2)

            //    //double meanI = 0;
            //    //double varI = 0;
            //    //for(int i = 0; i < corrlength; i++) {
            //    //    meanI += vals[i];
            //    //}
            //    //meanI /= vals.Length;
            //    //for(int i = 0; i < corrlength; i++) {
            //    //    varI += Math.Pow(vals[i] - meanI, 2);
            //    //}
            //    //varI /= vals.Length;
            //    //double[] autocorr = new double[corrlength];
            //    //for(int i = 0; i < vals.Length - corrlength; i++) {
            //    //    for(int j = 0; j < corrlength - i; j++) {
            //    //        autocorr[j] += (vals[i] - meanI) * (vals[i + j] - meanI);
            //    //    }
            //    //}
            //    //for(int i = 0; i < corrlength; i++) autocorr[i] /= varI;
            //    //rwplot.AddLine("Autocorrelation O(n²)", OxyColors.Black);
            //    //for(int i = 0; i < corrlength; i++) rwplot.AddPoint(i, autocorr[i] / (vals.Length - corrlength));




            //    //Calculate Autocorrelation with Wiener Kinchim Theorem O(n log n)

            //    Complex[] samples = new Complex[vals.Length * 2];
            //    for(int i = 0; i < vals.Length; i++) samples[i] = new Complex(vals[i], 0);
            //    for(int i = vals.Length; i < samples.Length; i++) samples[i] = new Complex();


            //    Transform.FourierForward(samples);
            //    //rwplot.AddLine("Wiener Khinchim FT", OxyColors.Automatic);
            //    //for(int i = 0; i < samples.Length / 2; i++) rwplot.AddPoint(i, samples[i].Real / 72);
            //    //rwplot.AddLine("Wiener Khinchim Imaginary FT", OxyColors.Automatic);
            //    //for(int i = 0; i < samples.Length / 2; i++) rwplot.AddPoint(i, samples[i].Imaginary / 72);


            //    for(int i = 0; i < samples.Length; i++) samples[i] = new Complex(samples[i].Magnitude * samples[i].Magnitude, 0);
            //    Transform.FourierInverse(samples);

            //    rwplot.AddLine("Wiener Khinchim", OxyColors.Automatic);
            //    for(int i = 0; i < samples.Length / 2; i++) rwplot.AddPoint(i, samples[i].Real / 72);
            //    rwplot.AddLine("Wiener Khinchim Imaginary", OxyColors.Automatic);
            //    for(int i = 0; i < samples.Length / 2; i++) rwplot.AddPoint(i, samples[i].Imaginary / 72);
            //}




            // 7

            //rwplot = new RandomWalkPlot(window.Plot1, "Kinetic Monte Carlo", "Fraction of particles in state 2 after different time");

            //{//7.1b
            //    int[] times = { 50, 300, 3000 };
            //    foreach(var t in times) {
            //        rwplot.AddLine("Analytic solution t = " + (t > 5000 ? "\u221E" : t.ToString()), OxyColors.Automatic);
            //        for(double eps = 0; eps <= 15; eps += 0.1) {
            //            double[] p = { 1, 0, 0 };
            //            rwplot.AddPoint(eps, timestep(p, eps, t));
            //        }
            //    }
            //}



            //{//real KMC 7.1c
            //    int N = 100000;
            //    int[] times = { 50, 300, 3000 };
            //    foreach(var t in times) {
            //        rwplot.AddLine("KMC t = " + (t > 5000 ? "\u221E" : t.ToString()), OxyColors.Automatic);
            //        for(double eps = 0; eps <= 15; eps += 0.2) {
            //            var erg = new double[3];
            //            for(int i = 0; i < N; i++) {
            //                if(i % 1000 == 0) { Console.WriteLine(eps + " " + t + " " + i / 1000 + "k"); }
            //                erg[kinmc(t, eps)]++;
            //            }
            //            rwplot.AddPoint(eps, (erg[2]) / N);
            //        }
            //    }
            //}

            //{//7.1d

            //    rwplot = new RandomWalkPlot(window.Plot1, "Kinetic Monte Carlo", "Cumulative bond traffic (T) and flux (K) for different energy values");

            //    int N = 5000;
            //    int[] epss = { 6, 7, 8, 10, 15 };
            //    foreach(var eps in epss) {
            //        long[,] kpm = new long[1000 / CL, 2];
            //        for(int i = 0; i < N; i++) {
            //            if(i % 1000 == 0) { Console.WriteLine(eps + " " + i / 10000 + "0k"); }
            //            kinmcflux(kpm, eps);
            //        }
            //        for(int i = 1; i < kpm.Length / 2; i++) {
            //            kpm[i, 0] += kpm[i - 1, 0];
            //            kpm[i, 1] += kpm[i - 1, 1];
            //        }

            //        bool TLOG = !true;
            //        byte cb = (byte)(100 + eps * 10);
            //        byte cbi = (byte)(255 - cb);
            //        int x = 0;
            //        if(x != 2) {
            //            rwplot.AddLine("ln(T), ε  = " + eps, OxyColor.FromArgb(cb, 0, 0, cb));
            //            for(int i = 0; i < kpm.Length / 2; i++) {
            //                double T = kpm[i, 0] + kpm[i, 1];
            //                double y = Math.Log(i * CL);
            //                if(!TLOG) y = i * CL;
            //                rwplot.AddPoint(y, Math.Log(T / N));
            //            }
            //        }
            //        if(x != 1) {
            //            rwplot.AddLine("K, ε  = " + eps, OxyColor.FromArgb(cb, 0, cb, cb));
            //            for(int i = 0; i < kpm.Length / 2; i++) {
            //                double k = kpm[i, 0] - kpm[i, 1];
            //                double y = Math.Log(i * CL);
            //                if(!TLOG) y = i * CL;
            //                rwplot.AddPoint(y, k / N);
            //            }
            //        }

            //    }
            //}

            {//7.2a
                int[] css = { 2,5,10 };
                foreach(var c in css)
                    //geneticSwitch(c, 1, 1, 50000, true, false);
                ;
                int[] cs = { 2, 5, 10 };
                foreach(var c in cs)
                    geneticSwitch(c, 0.5, 1, 1, false, false);
                    ;

            }
        }
        int CL = 2; // 7.1


        private int kinmc(int endtime, double eps) {
            double m = 10;
            double c = 0.01;
            double time = 0;
            int state = 0, nstate = 0;
            double r01 = c * m;
            double r02 = c;
            double r1 = Math.Exp(-eps / 2);
            double r2 = Math.Exp(-eps);

            double[] rtot = { r01 + r02, r1, r2 };

            double[] swtime = new double[3];
            while(time <= endtime) {
                state = nstate;
                if(state == 0)
                    nstate = rng.NextDouble() * rtot[state] > r01 ? 2 : 1;
                else
                    nstate = 0;
                time += expverteilung(rtot[state]);
            }
            return state;
        }

        private void kinmcflux(long[,] kpm, double eps) {
            double m = 10;
            double c = 0.01;

            double r01 = c * m;
            double r02 = c;
            double r1 = Math.Exp(-eps / 2);
            double r2 = Math.Exp(-eps);

            double endtime = 1000;
            double time = 0;


            int state = 0, nstate = 0;

            double[] rtot = { r01 + r02, r1, r2 };
            double[] swtime = new double[3];
            while(time <= endtime) {
                state = nstate;
                int itime = ((int)time) / CL;
                switch(state) {
                    case 0: nstate = rng.NextDouble() * rtot[state] > r01 ? 2 : 1;
                        kpm[itime, 0] += nstate;
                        break;
                    default: nstate = 0;
                        kpm[itime, 1] += state;
                        break;
                }
                time += expverteilung(rtot[state]);
            }
            //if(state == 0)
            //    kpm[kpm.Length/2 - 1, 0] -= nstate;
            //else
            //    kpm[kpm.Length/2 - 1, 1] += state;

        }


        enum oState { e, A, B };
        enum oAction { pA, pB, cA, cB, uC, mA, mB, pA2, pB2, mA2, mB2 }
        oAction[] actions = { oAction.pA, oAction.pB, oAction.cA, oAction.cB, oAction.uC, oAction.mA, oAction.mB, oAction.pA2, oAction.pB2, oAction.mA2, oAction.mB2 };
        private void geneticSwitch(int c, double k, int N, double endtime, bool detailedpath, bool verydetail) {
            rwplot = rwplot ?? new RandomWalkPlot(window.Plot1, "Kinetic Monte Carlo", "Genetic switch k=" + k + "; c=" + c + "; " + N + " trajector" + (N > 1 ? "ies" : "y") + "; maxTime = " + endtime);
            int[] nA = new int[N], nB = new int[N], nA2 = new int[N], nB2 = new int[N];
            bool[] cont = Enumerable.Repeat<bool>(true, N).ToArray();
            var time = new double[N];
            oState[] state = new oState[N];
            rwplot.SetAxisTitles("λ", "p(λ)");


            LineSeries pnA = null;
            LineSeries pnB = null;
            LineSeries pnA2 = null;
            LineSeries pnB2 = null;
            LineSeries pL = null;
            if(detailedpath) {
                rwplot.SetAxisTitles("Time", "Amount");
                if(verydetail) {

                    pnA = rwplot.AddLineM("A, c = " + c, OxyColors.Automatic);
                    pnB = rwplot.AddLineM("B, c = " + c, OxyColors.Automatic);
                    pnA2 = rwplot.AddLineM("A2, c = " + c, OxyColors.Automatic);
                    pnB2 = rwplot.AddLineM("B2, c = " + c, OxyColors.Automatic);
                } else
                    pL = rwplot.AddLineM("λ, c = " + c, OxyColors.Automatic);
            }
            double[] lambda = new double[N];
            int finished = 0;
            while(finished != N) {
                for(int s = 0; s < N; s++) {
                    if(cont[s]) {
                        var rates = getRates(nA[s], nA2[s], nB[s], nB2[s], state[s], c, k);
                        double rand = rng.NextDouble() * rates[11];
                        int actindex = 0;
                        double comp = rates[0];
                        while(rand > comp) {
                            actindex++;
                            comp += rates[actindex];
                        }
                        double dt = expverteilung(rates[11]);

                        oAction todo = actions[actindex];
                        if(state[s] == oState.e && todo == oAction.uC)
                            Console.WriteLine("LOLWTF");
                        if(state[s] == oState.A && todo == oAction.pB)
                            Console.WriteLine("LOLWTF");
                        if(state[s] == oState.B && todo == oAction.pA)
                            Console.WriteLine("LOLWTF");
                        if(state[s] != oState.e && (todo == oAction.cB || todo == oAction.cA))
                            Console.WriteLine("LOLWTF");
                        if(nA[s] == 0 && (todo == oAction.mA || todo == oAction.pA2))
                            Console.WriteLine("LOLWTF");
                        if(nB[s] == 0 && (todo == oAction.mB || todo == oAction.pB2))
                            Console.WriteLine("LOLWTF");
                        if(nA2[s] == 0 && (todo == oAction.mA2 || todo == oAction.cA))
                            Console.WriteLine("LOLWTF");
                        if(nB2[s] == 0 && (todo == oAction.mB2 || todo == oAction.cB))
                            Console.WriteLine("LOLWTF");

                        time[s] += dt;
                        if(time[s] > endtime && cont[s]) {
                            cont[s] = false;
                            finished++;
                            continue;
                        }
                        switch(actions[actindex]) {
                            case oAction.pA: nA[s]++; break;
                            case oAction.pB: nB[s]++; break;
                            case oAction.cA: state[s] = oState.A; nA2[s] -= 1; break;
                            case oAction.cB: state[s] = oState.B; nB2[s] -= 1; break;
                            case oAction.uC:
                                if(state[s] == oState.A)
                                    nA2[s] += 1;
                                else
                                    nB2[s] += 1;
                                state[s] = oState.e;
                                break;
                            case oAction.mA: nA[s]--; break;
                            case oAction.mB: nB[s]--; break;
                            case oAction.pA2: nA2[s]++; nA[s] -= 2; break;
                            case oAction.pB2: nB2[s]++; nB[s] -= 2; break;
                            case oAction.mA2: nA2[s]--; nA[s] += 2; break;
                            case oAction.mB2: nB2[s]--; nB[s] += 2; break;
                        }
                        lambda[s] = (nA[s] + 2 * nA2[s]) - (nB[s] + 2 * nB2[s]);
                        if(detailedpath) {
                            if(verydetail) {
                                pnA.Points.Add(new DataPoint(time[s], nA[s]));
                                pnB.Points.Add(new DataPoint(time[s], nB[s]));
                                pnA2.Points.Add(new DataPoint(time[s], nA2[s]));
                                pnB2.Points.Add(new DataPoint(time[s], nB2[s]));
                            } else
                                pL.Points.Add(new DataPoint(time[s], lambda[s]));
                        }
                        if(nA[s] < 0 || nB[s] < 0)
                            Console.WriteLine("LOLWTF");
                    }
                }
            }
            if(!detailedpath) {
                LineSeries pLambda = rwplot.AddLineM("λ, c = " + c, OxyColors.Automatic);
                double lp = lambda.Max();
                double lm = lambda.Min();
                int classes = 50;
                double di = 2 * lp / classes;
                di = 1;
                for(double i = lm; i < lp; i += di) {
                    double n = 0;
                    foreach(var y in lambda) {
                        if(i == y)
                            n++;
                    }
                    pLambda.Points.Add(new DataPoint(i, n / N));
                }
            }
        }

        private double[] getRates(int na, int na2, int nb, int nb2, oState state, int c, double k) {
            double[] rates = { 
                                 k, //pA 0
                                 k, //pB 1
                                 k*c*na2, //cA 2
                                 k*c*nb2, //cB 3
                                 0, //uC 4
                                 k*((double)na)*0.25, //mA 5 
                                 k*((double)nb)*0.25, //mB 6 
                                 k*(na/2), //pA2 7
                                 k*(nb/2), //pB2 8 
                                 k*((double)na2)*5, //mA2 9
                                 k*((double)nb2)*5, //mB2 10
                                 0 // total rate - calculated later 11
                             };
            switch(state) {
                case oState.A:
                    rates[1] = 0;
                    rates[2] = 0;
                    rates[3] = 0;
                    rates[4] = k;
                    break;
                case oState.B:
                    rates[0] = 0;
                    rates[2] = 0;
                    rates[3] = 0;
                    rates[4] = k;
                    break;
            }
            if(na == 1)
                rates[7] = 0;
            if(nb == 1)
                rates[8] = 0;
            rates[11] = rates.Sum();
            return rates;
        }


        private double timestep(double[] p, double eps, int t) {
            double c = 0.01;
            double m = 10;
            for(int i = 0; i < t; i++) {
                p[0] = p[0] * (1 - c - c * m) + p[1] * Math.Exp(-eps) + p[2] * Math.Exp(-eps / 2);
                p[1] = p[0] * c + p[1] * (1 - Math.Exp(-eps));
                p[2] = 1 - p[0] - p[1];
            }
            return p[1];
        }


        private double minX = -10;
        private double scale = 0.1;
        private void TakeWalk(Walker walker, OxyColor col) {
            for(int traj = 0; traj < Traj; traj++) {
                if(traj % 1000 == 0)
                    Console.WriteLine(traj);
                walker.Reset();
                for(int step = 1; step <= 10; step++) {
                    double pos = walker.Walk(step * 2 - 1);
                    for(int i = 0; i < classesHist; i++) {
                        double classL = minX + i * scale;
                        if(pos >= classL && pos <= classL + scale) {
                            positions[i, step]++;
                            break;
                        }
                    }
                }
            }
        }
        //var s21 = new RectangleBarSeries();
        //var s10 = new RectangleBarSeries();
        //s10.Title = "Histogram Cauchy Distribution (scaled to fit)";
        ////s21.Title = "Histogram n = 10";

        //double mean = 0;
        //for(int i = 0; i < Traj; i++) {
        //    position[i] /= Steps;
        //    mean += position[i];
        //}
        //mean /= Traj;
        //Console.WriteLine(mean);
        //int classN = 100;
        //double intervall = 10;
        //double isteps = 2 * intervall / classN;
        //double[] classes = new double[classN];
        //for(int i = 0; i < Traj; i++) {
        //    //Console.WriteLine(position[i]);
        //    for(int c = 0; c < classN; c++) {
        //        double classMin = -intervall + c * isteps;
        //        if(position[i] >= classMin && position[i] < classMin + isteps) {
        //            classes[c] += 5;
        //        }
        //    }
        //}
        //for(int c = 0; c < classN; c++) {
        //    double classMin = -intervall + c * isteps;
        //    s10.Items.Add(new RectangleBarItem() { X0 = classMin, X1 = classMin + isteps, Y0 = 0, Y1 = classes[c] });
        //}

        //window.Plot1.Series.Add(s10);

        private double heatEq(double x, double t) {
            double a = 0.5;
            return Math.Exp(-(x * x) / (4 * a * t)) / Math.Sqrt(Math.PI * a * t);
        }
        private double cauchyverteilung(double x) {
            return 1 / (Math.PI * (1 + x * x));
        }
        private double normalverteilung(double x, double mu, double varianz) {
            return Math.Exp(-Math.Pow(x - mu, 2) / (2 * varianz)) / Math.Sqrt(2 * Math.PI * varianz);
        }
        private double reflektierteverteilung(double D, double x, double t, double lowerbound, double upperbound) {
            double L = upperbound - lowerbound;
            double sum = 0;
            for(int n = 1; n < 2000; n++) {
                sum += Math.Cos(2 * Math.PI * n * x / L) * Math.Exp(-4 * Math.PI * Math.PI * n * n * D * t / L / L);
            }
            return 1 / L + 2 / L * sum;
        }
        private double reflektiertmsd(double D, double t, double L) {
            //return 1.0 / 12 - 1 / Math.PI / Math.PI * Math.Exp(-4 * Math.PI * Math.PI * D * t / L / L);
            double sum = 0;
            for(int n = 1; n < 5000; n++) {
                sum += Math.Pow(-1, n) / n / n * Math.Exp(-4 * Math.PI * Math.PI * n * n * D * t / L / L);
            }
            return (1.0 / 12) + (sum / Math.PI / Math.PI);
        }

        private double expverteilung(double rate) {
            return -1 / rate * Math.Log(rng.NextDouble());
        }
        private readonly double sqrt3 = Math.Sqrt(3);
        private double SumUniform(int n) {
            double y = 0;
            for(int i = 0; i < n; i++) {
                double x = rng.NextDouble() * sqrt3;
                if(rng.NextDouble() < 0.5)
                    y -= x;
                else
                    y += x;
            }
            return y / Math.Sqrt(n);
        }

        private int acc = 0;
        private int rej = 0;
        double variance = 0.295;
        double multi = 1.4;
        private double DellagoRandom() {
            acc++;
        start:
            rej++;
            double x = getGaussian() * Math.Sqrt(variance);
            if(rng.NextDouble() < 0.5)
                x -= 1;
            else
                x += 1;
            if(rng.NextDouble() * DellagoMajor(x) < DellagoDistrib(x)) return x;
            goto start;

        }
        private double DellagoDistrib(double x) {
            return 0.50665436034 * Math.Exp(-(x * x - 1) * (x * x - 1));
        }

        private double DellagoMajor(double x) {
            return 0.5 * multi * (normalverteilung(x, -1, variance) + normalverteilung(x, 1, variance));
        }
        private double getGaussian() {
            double u1 = rng.NextDouble();
            double u2 = rng.NextDouble();
            return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2.0 * Math.PI * u2);
        }
    }
    public class HistogramBin {
        public string Label { get; set; }
        public double Value { get; set; }
    }
}
