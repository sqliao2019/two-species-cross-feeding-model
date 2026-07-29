input.Sin = 5;
input.Min = 0.0;

Param = makeSym(getPara(5));

ScanParam1D(Param, 'k', [0.1 0.4], 10000);

ScanParam1D(Param, 'dilution', [0.0 2.0], 10000);

Param.dilution = 1.5;
ScanParam1D(Param, 'Min', [0.0 0.8], 10000);
