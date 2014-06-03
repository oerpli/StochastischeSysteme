using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using OxyPlot;
using OxyPlot.Axes;
using OxyPlot.Series;
using StochSyst.Annotations;

namespace StochSyst
{
    public class BoxPlot
    {
        public List<LineSeries> lines { private get; set; }
        private LineSeries lastLine = null;
        private PlotModel pm;
        public BoxPlot(PlotModel pm, String title, String subTitle)
        {
            this.pm = pm;
            pm.LegendSymbolLength = 10;
            pm.LegendPosition = LegendPosition.TopLeft;
            pm.Title = title;
            pm.Subtitle = subTitle;
            var linearAxis1 = new LinearAxis();
            var linearAxis2 = new LinearAxis();
            linearAxis1.Position = AxisPosition.Bottom;
            pm.Axes.Add(linearAxis1);
            pm.Axes.Add(linearAxis2);
        }

        public BoxPlot AddLine(String lineTitle, OxyColor? color)
        {
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

        public BoxPlot AddPoint(double x, double y)
        {
            (lastLine ?? AddLine("Untitled", lastLine.Color).lastLine).Points.Add(new DataPoint(x, y));
            return this;
        }
    }
}
