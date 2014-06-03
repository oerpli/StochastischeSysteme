using System;
using OxyPlot.Series;

public abstract class Walker {
    static protected Random rng = new Random();
    protected double pos = 0;
    public double Walk(int steps) {
        for(int i = 0; i < steps; i++) {
            NextStep();
        }
        return pos;
    }
    protected abstract void NextStep();
    public abstract string Name();
    internal void Reset() {
        pos = 0;
    }
}

public class WalkerCauchy : Walker {
    protected override void NextStep() {
        pos += cauchyRandom();
    }
    public override string Name() {
        return "Cauchy";
    }
    private double cauchyRandom() {
        return (Math.Tan((rng.NextDouble() - 0.5) * Math.PI));
    }
}
public class WalkerBinomial : Walker {
    public override string Name() {
        return "Binomial";
    }
    protected override void NextStep() {
        if(rng.NextDouble() < 0.5)
            pos++;
        else
            pos--;
    }
}
public class WalkerGaussian : Walker {
    public override string Name() {
        return "Gaussian";
    }
    protected override void NextStep() {//using Box-Muller transform 
        pos += getGaussian();
    }
    protected double getGaussian() {
        double u1 = rng.NextDouble();
        double u2 = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2.0 * Math.PI * u2);
    }
}

public class WalkerGaussianReflective : WalkerGaussian {
    private byte handling = 0;
    private double old = 0;
    private double LU = 10;
    private double LO = -10;
    public WalkerGaussianReflective(byte handling, double lowerBound, double upperBound)
        : base() {
        this.handling = handling;
        this.LO = lowerBound;
        this.LU = upperBound;
    }
    protected override void NextStep() {
        old = pos;
        pos += getGaussian();
        if(pos <= LO || pos >= LU) {
            switch(handling) {
                case 0: ReflectStep(); return;
                case 1: TruncateStep(); return;
                case 2: DismissStep(); return;
                case 3: RetryStep(); return;
            }
        }
    }
    private void ReflectStep() {
        if(pos <= LO)
            pos = pos - 2 * (pos - LO);
        if(pos >= LU)
            pos = pos - 2 * (pos - LU);
    }
    private void TruncateStep() {
        pos = Math.Max(LO, Math.Min(pos, LU));
    }
    private void DismissStep() {
        pos = old;
    }
    private void RetryStep() {
        while(pos < LO || pos > LU) {
            pos = old + getGaussian();
        }
    }
}