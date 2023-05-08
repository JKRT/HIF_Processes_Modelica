package ThetaPaperExample

model M
  Real x1;
  Real x2;
  Real x3;
  Real x4;
  Real x5;
  Real u;
equation
  u = time;
  der(x1) = u;
  der(x2) = x1 + x2;
  der(x3) = x2 * x4 + 1;
  der(x4) = x3 + 2 * x5 + 2;
  der(x5) = x4 + 3 * x5 + 3;
end M;

end ThetaPaperExample;  
