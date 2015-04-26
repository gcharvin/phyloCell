function P=computeProba(pdfout,range,enable,varz)


%% calculate effective proba knowing kernel density estimators
% input X,Y,A,I,DY,DA,DI

pDX = evaluatePointsUnderPdf(pdfout(1), [varz(1) ;varz(3); varz(5)]);
pDY = evaluatePointsUnderPdf(pdfout(2), [varz(2) ;varz(3); varz(6)]);
pDA = evaluatePointsUnderPdf(pdfout(3), [varz(4) ;varz(3); varz(7)]);
pDI = evaluatePointsUnderPdf(pdfout(4), [varz(4) ;varz(3); varz(8)]);


P=1;

if enable(1)
    P=P*pDX*range(1)*range(3);
end
if enable(2)
    P=P*pDY*range(2)*range(3);
end
if enable(3)
    P=P*pDA*range(3)*range(4);
end
if enable(4)
    P=P*pDI*range(4)*range(3);
end

