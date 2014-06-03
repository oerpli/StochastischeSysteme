using System;
class Quine {
    static string f =
        "using System;class Quine{{static string f={1}{0}{1};static void Main(){{Console.Write(f,f,Convert.ToChar(34));}}}}";
    static void Main() { Console.Write(f, f, Convert.ToChar(34)); }
}