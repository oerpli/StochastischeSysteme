using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using StochSyst.ViewModels;
using System.Globalization;


using OxyPlot;
using OxyPlot.Series;


namespace StochSyst
{
    abstract class Prob1
    {
        public static void main(MainWindowModel windowModel)
        {
            new RandomWalk(windowModel);
            //new Sampler(windowModel);
            //new RandomWalk2D(windowModel, "weird_walker.txt", "walker in 2d space", 1, 800000);
            //new RandomWalk2D(windowModel, "strange_walkers.txt", "walkers in 2d space", 1, 2500000);
        }
    }
}
