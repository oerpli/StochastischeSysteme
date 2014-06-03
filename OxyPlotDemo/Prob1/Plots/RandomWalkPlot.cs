using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using OxyPlot;
using OxyPlot.Axes;
using OxyPlot.Series;
using StochSyst.Annotations;

namespace StochSyst {
    public class RandomWalkPlot {
        public List<LineSeries> lines { private get; set; }
        private LineSeries lastLine = null;
        private PlotModel pm;
        private LinearAxis xaxe;
        private LinearAxis yaxe;
        public RandomWalkPlot(PlotModel pm, String title, String subTitle) {
            this.pm = pm;
            pm.LegendSymbolLength = 30;
            pm.LegendPosition = LegendPosition.TopLeft;
            pm.Title = title;
            pm.Subtitle = subTitle;
            xaxe = new LinearAxis();
            yaxe = new LinearAxis();
            xaxe.Position = AxisPosition.Bottom;
            pm.Axes.Add(xaxe);
            pm.Axes.Add(yaxe);
        }
        public void SetAxisTitles(String x, String y) {
            xaxe.Title = x;
            yaxe.Title = y;
        }

        public RandomWalkPlot AddLine(String lineTitle, OxyColor? color) {
            lastLine = new LineSeries();
            lastLine.Color = color ?? OxyColors.Black;
            //lastLine.MarkerFill = OxyColors.MidnightBlue;
            //lastLine.MarkerSize = 2;
            //lastLine.MarkerStroke = OxyColors.Navy;
            //lastLine.MarkerStrokeThickness = 0.25;
            //lastLine.MarkerType = MarkerType.Circle;
            lastLine.CanTrackerInterpolatePoints = false;
            lastLine.Title = lineTitle;
            pm.Series.Add(lastLine);
            return this;
        }
        public LineSeries AddLineM(String lineTitle, OxyColor? color) {
            lastLine = new LineSeries();
            lastLine.Color = color ?? OxyColors.Black;
            //lastLine.MarkerFill = OxyColors.MidnightBlue;
            //lastLine.MarkerSize = 2;
            //lastLine.MarkerStroke = OxyColors.Navy;
            //lastLine.MarkerStrokeThickness = 0.25;
            //lastLine.MarkerType = MarkerType.Circle;
            lastLine.CanTrackerInterpolatePoints = false;
            lastLine.Title = lineTitle;
            pm.Series.Add(lastLine);
            return lastLine;
        }

        public RandomWalkPlot AddPoint(double x, double y) {
            (lastLine ?? AddLine("Untitled", lastLine.Color).lastLine).Points.Add(new DataPoint(x, y));
            return this;
        }
        public RandomWalkPlot AddPoint(LineSeries tar,double x, double y) {
            (tar ?? AddLine("Untitled", lastLine.Color).lastLine).Points.Add(new DataPoint(x, y));
            return this;
        }
    }
}
