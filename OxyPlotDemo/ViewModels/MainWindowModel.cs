using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using OxyPlot;
using OxyPlot.Axes;
using OxyPlot.Series;
using StochSyst.Annotations;

namespace StochSyst.ViewModels
{
    public class MainWindowModel : INotifyPropertyChanged
    {
        private PlotModel plot1;
        public PlotModel Plot1 { get { return plot1; } set { plot1 = value; OnPropertyChanged("Plot1"); } }
        //  public PlotModel Plot2 { get { return plot2; } set { plot2 = value; OnPropertyChanged("Plot2"); } }
        public String PlotTitle = "Plot Title!";

        private DateTime lastUpdate = DateTime.Now;

        public MainWindowModel()
        {
            Plot1 = new PlotModel();
            Prob1.main(this);
            //SetUpModel(Plot1);
            //SetUpModel(Plot2);
            //LoadData(Plot1);
            //LoadData(Plot2);
            //Plot1.Series.Add(new FunctionSeries(Math.Cos, 0, 100, 0.1, "cos(x)"));
        }

        private readonly List<OxyColor> colors = new List<OxyColor> { OxyColors.AliceBlue, OxyColors.IndianRed, OxyColors.Coral, OxyColors.Chartreuse, OxyColors.Azure };
        private readonly List<MarkerType> markerTypes = new List<MarkerType> { MarkerType.Plus, MarkerType.Star, MarkerType.Diamond, MarkerType.Triangle, MarkerType.Cross };

        private void SetUpModel(PlotModel pm)
        {
            pm.LegendTitle = "Legend";
            pm.LegendOrientation = LegendOrientation.Horizontal;
            pm.LegendPlacement = LegendPlacement.Outside;
            pm.LegendPosition = LegendPosition.TopRight;
            pm.LegendBackground = OxyColor.FromAColor(200, OxyColors.White);
            pm.LegendBorder = OxyColors.Black;

            var dateAxis = new DateTimeAxis(AxisPosition.Bottom, "Date", "HH:mm") { MajorGridlineStyle = LineStyle.Solid, MinorGridlineStyle = LineStyle.Dot, IntervalLength = 80 };
            pm.Axes.Add(dateAxis);
            var valueAxis = new LinearAxis(AxisPosition.Left, 0) { MajorGridlineStyle = LineStyle.Solid, MinorGridlineStyle = LineStyle.Dot, Title = "Value" };
            pm.Axes.Add(valueAxis);

        }

        private void LoadData(PlotModel pm)
        {
            List<Measurement> measurements = Data.GetData();
            var dataPerDetector = measurements.GroupBy(m => m.DetectorId).OrderBy(m => m.Key).ToList();
            foreach (var data in dataPerDetector)
            {
                var lineSerie = new LineSeries
                {
                    StrokeThickness = 2,
                    MarkerSize = 3,
                    MarkerStroke = colors[data.Key],
                    MarkerType = markerTypes[data.Key],
                    CanTrackerInterpolatePoints = false,
                    Title = string.Format("Detector {0}", data.Key),
                    Smooth = false,
                };

                data.ToList().ForEach(d => lineSerie.Points.Add(new DataPoint(DateTimeAxis.ToDouble(d.DateTime), d.Value)));
                pm.Series.Add(lineSerie);
            }
            lastUpdate = DateTime.Now;
        }

        public void UpdateModel(PlotModel pm)
        {
            List<Measurement> measurements = Data.GetUpdateData(lastUpdate);
            var dataPerDetector = measurements.GroupBy(m => m.DetectorId).OrderBy(m => m.Key).ToList();
            foreach (var data in dataPerDetector)
            {
                var lineSerie = Plot1.Series[data.Key] as LineSeries;
                if (lineSerie != null)
                {
                    data.ToList()
                        .ForEach(d => lineSerie.Points.Add(new DataPoint(DateTimeAxis.ToDouble(d.DateTime), d.Value)));
                }
            }
            lastUpdate = DateTime.Now;
        }

        public event PropertyChangedEventHandler PropertyChanged;

        [NotifyPropertyChangedInvocator]
        protected virtual void OnPropertyChanged(string propertyName)
        {
            PropertyChangedEventHandler handler = PropertyChanged;
            if (handler != null) handler(this, new PropertyChangedEventArgs(propertyName));
        }

        internal void UpdateModel()
        {
            UpdateModel(Plot1);
        }
    }
}
