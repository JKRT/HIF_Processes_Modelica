package TestThetaMethod "Test cases for Dirk Zimmer's Theta Method"

model ThetaCapacitator
  extends Modelica.Electrical.Analog.Interfaces.OnePort(v(start=0));
  parameter Modelica.SIunits.Capacitance C(start=1) "Capacitance";
  parameter Real THETA;
equation
  i = C*THETA*der(v);
annotation(
      __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian -d=infoXmlOperations ",
      __OpenModelica_simulationFlags(lv = "LOG_STATS", s = "dassl"),
  Diagram(graphics = {Line(origin = {0, 58.1667}, points = {{0, -39.5}, {0, 40.5}, {0, 38.5}}), Line(origin = {0, 58.1667}, points = {{0, -39.5}, {0, 40.5}, {0, 38.5}}), Line(origin = {0, -59}, points = {{0, -39}, {0, 39}})}));
end ThetaCapacitator;

  extends Modelica.Icons.Package;

  package NonlinearCircuit
    extends Modelica.Icons.Package;

    package Test

      extends Modelica.Icons.ExamplesPackage;

      model Circuit1Static
        extends Modelica.Icons.Example;
        Modelica.Electrical.Analog.Basic.Ground ground annotation(
          Placement(visible = true, transformation(origin = {0, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Electrical.Analog.Sources.StepCurrent stepCurrent(I = 1) annotation(
          Placement(visible = true, transformation(origin = {-80, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
        Modelica.Electrical.Analog.Basic.Capacitor C1(C(displayUnit = "uF") = 0.0001000000000000001) annotation(
          Placement(visible = true, transformation(origin = {-54, -10}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
        Modelica.Electrical.Analog.Basic.Resistor R1(R = 1000) annotation(
          Placement(visible = true, transformation(origin = {78, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Electrical.Analog.Semiconductors.Diode D1(Ids = 1e-9, Maxexp = 40)  annotation(
          Placement(visible = true, transformation(origin = {-12, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Electrical.Analog.Semiconductors.Diode D2(Ids = 1e-9, Maxexp = 40)  annotation(
          Placement(visible = true, transformation(origin = {-12, 24}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Electrical.Analog.Semiconductors.Diode D3(Ids = 1e-9, Maxexp = 40) annotation(
          Placement(visible = true, transformation(origin = {22, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Electrical.Analog.Semiconductors.Diode diode(Ids = 1e-9, Maxexp = 40) annotation(
          Placement(visible = true, transformation(origin = {50, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(C1.n, ground.p) annotation(
          Line(points = {{-54, -20}, {-54, -30}, {0, -30}}, color = {0, 0, 255}));
        connect(stepCurrent.p, ground.p) annotation(
          Line(points = {{-80, -20}, {-80, -30}, {0, -30}}, color = {0, 0, 255}));
        connect(R1.n, ground.p) annotation(
          Line(points = {{88, 10}, {96, 10}, {96, -30}, {0, -30}}, color = {0, 0, 255}));
        connect(stepCurrent.n, C1.p) annotation(
          Line(points = {{-80, 0}, {-80, 10}, {-54, 10}, {-54, 0}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
        connect(C1.p, D1.p) annotation(
          Line(points = {{-54, 0}, {-54, 10}, {-32, 10}, {-32, 0}, {-22, 0}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
        connect(C1.p, D2.n) annotation(
          Line(points = {{-54, 0}, {-54, 10}, {-32, 10}, {-32, 24}, {-22, 24}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
        connect(D2.p, D3.p) annotation(
          Line(points = {{-2, 24}, {6, 24}, {6, 10}, {12, 10}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
        connect(D1.n, D3.p) annotation(
          Line(points = {{-2, 0}, {6, 0}, {6, 10}, {12, 10}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
        connect(D3.n, diode.p) annotation(
          Line(points = {{32, 10}, {40, 10}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
        connect(diode.n, R1.p) annotation(
          Line(points = {{60, 10}, {68, 10}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
        annotation(
          Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}})),
  __OpenModelica_simulationFlags(lv = "LOG_SIMULATION,LOG_STATS", noEquidistantTimeGrid = "()", s = "dassl", cpu = "()"),
  __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian -d=infoXmlOperations ");
      end Circuit1Static;
      model Circuit1Dynamic
        extends Circuit1Static;
  Modelica.Electrical.Analog.Basic.Capacitor Cnl(C (displayUnit = "F")= 1e-12) annotation(
          Placement(visible = true, transformation(origin = {60, -10}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
        equation
          connect(Cnl.n, ground.p) annotation(
          Line(points = {{60, -20}, {60, -30}, {0, -30}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
          connect(diode.n, Cnl.p) annotation(
          Line(points = {{60, 10}, {60, 0}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
        annotation(
          Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}})),
  __OpenModelica_simulationFlags(lv = "LOG_STATS", s = "dassl", cpu = "()"),
  __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian");
      end Circuit1Dynamic;

model ThetaCircuit1Static
         extends Modelica.Icons.Example;
        Modelica.Electrical.Analog.Basic.Ground ground annotation(
          Placement(visible = true, transformation(origin = {0, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Electrical.Analog.Sources.StepCurrent stepCurrent(I = 1) annotation(
          Placement(visible = true, transformation(origin = {-80, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
        ThetaCapacitator C1(C(displayUnit = "uF") = 0.0001000000000000001) annotation(
          Placement(visible = true, transformation(origin = {-54, -10}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
        Modelica.Electrical.Analog.Basic.Resistor R1(R = 1000) annotation(
          Placement(visible = true, transformation(origin = {78, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Electrical.Analog.Semiconductors.Diode D1(Ids = 1e-9, Maxexp = 40)  annotation(
          Placement(visible = true, transformation(origin = {-12, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Electrical.Analog.Semiconductors.Diode D2(Ids = 1e-9, Maxexp = 40)  annotation(
          Placement(visible = true, transformation(origin = {-12, 24}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Electrical.Analog.Semiconductors.Diode D3(Ids = 1e-9, Maxexp = 40) annotation(
          Placement(visible = true, transformation(origin = {22, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Electrical.Analog.Semiconductors.Diode diode(Ids = 1e-9, Maxexp = 40) annotation(
          Placement(visible = true, transformation(origin = {50, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(C1.n, ground.p) annotation(
    Line(points = {{-54, -20}, {-54, -30}, {0, -30}}, color = {0, 0, 255}));
  connect(stepCurrent.p, ground.p) annotation(
    Line(points = {{-80, -20}, {-80, -30}, {0, -30}}, color = {0, 0, 255}));
  connect(R1.n, ground.p) annotation(
    Line(points = {{88, 10}, {96, 10}, {96, -30}, {0, -30}}, color = {0, 0, 255}));
  connect(stepCurrent.n, C1.p) annotation(
    Line(points = {{-80, 0}, {-80, 10}, {-54, 10}, {-54, 0}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
  connect(C1.p, D1.p) annotation(
    Line(points = {{-54, 0}, {-54, 10}, {-32, 10}, {-32, 0}, {-22, 0}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
  connect(C1.p, D2.n) annotation(
    Line(points = {{-54, 0}, {-54, 10}, {-32, 10}, {-32, 24}, {-22, 24}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
  connect(D2.p, D3.p) annotation(
    Line(points = {{-2, 24}, {6, 24}, {6, 10}, {12, 10}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
  connect(D1.n, D3.p) annotation(
    Line(points = {{-2, 0}, {6, 0}, {6, 10}, {12, 10}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
  connect(D3.n, diode.p) annotation(
    Line(points = {{32, 10}, {40, 10}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
  connect(diode.n, R1.p) annotation(
    Line(points = {{60, 10}, {68, 10}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
  annotation(
    Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}})),
__OpenModelica_simulationFlags(lv = "LOG_SIMULATION,LOG_STATS", noEquidistantTimeGrid = "()", s = "dassl", cpu = "()"),
__OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian -d=infoXmlOperations -d=infoXmlOperations");
end ThetaCircuit1Static;

      model ThetaCircuit1Dynamic
        extends Circuit1Static;
      TestThetaMethod.ThetaCapacitator Cnl(C (displayUnit = "F")= 1e-12) annotation(
          Placement(visible = true, transformation(origin = {60, -14}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
        equation
          connect(Cnl.n, ground.p) annotation(
          Line(points = {{60, -24}, {60, -30}, {0, -30}}, color = {0, 0, 255}));
          connect(diode.n, Cnl.p) annotation(
          Line(points = {{60, 10}, {60, -4}}, color = {0, 0, 255}));
        annotation(
          Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}})),
      __OpenModelica_simulationFlags(lv = "LOG_DASSL,LOG_SIMULATION,LOG_STATS", s = "dassl", cpu = "()"),
      __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian -d=infoXmlOperations ");
      end ThetaCircuit1Dynamic;

      model ThetaCircuit2Dynamic
        parameter Real THETA = 1.0;
        extends Circuit1Static;
       Modelica.Electrical.Analog.Basic.Capacitor Cp(C (displayUnit = "F")= 1e-12 * THETA) annotation(
          Placement(visible = true, transformation(origin = {60, -14}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
        equation
          connect(Cp.n, ground.p) annotation(
          Line(points = {{60, -24}, {60, -30}, {0, -30}}, color = {0, 0, 255}));
          connect(diode.n, Cp.p) annotation(
          Line(points = {{60, 10}, {60, -4}}, color = {0, 0, 255}));
        annotation(
          Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}})),
      __OpenModelica_simulationFlags(lv = "LOG_DASSL,LOG_SIMULATION,LOG_STATS", s = "dassl", cpu = "()"),
      __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian -d=infoXmlOperations ");
      end ThetaCircuit2Dynamic;


    end Test;
  end NonlinearCircuit;
  annotation(
    uses(Modelica(version = "3.2.3")));
end TestThetaMethod;
