function [fitFunction, fitString] = fitArray(n,quantityArray,quantityEnd,fitType,startPoints)
% quantity end is rhoE, vE, uE etc
% fitType can be fourier8, sin8, general logistic, general logistic only
% upper bound

switch fitType
    case 'fourier8'
        % uses MATLAB curve fitting toolbox; 8 term fourier series works well at
        % low Mach, but at higher Mach less so; a different fit may be desirable
        fitFunctionBL = fit(n,quantityArray,'fourier8');

        % obtains coefficient values
        fitCoeffs = coeffvalues(fitFunctionBL);
        a0 = fitCoeffs(1);
        a1 = fitCoeffs(2);
        b1 = fitCoeffs(3);
        a2 = fitCoeffs(4);
        b2 = fitCoeffs(5);
        a3 = fitCoeffs(6);
        b3 = fitCoeffs(7);
        a4 = fitCoeffs(8);
        b4 = fitCoeffs(9);
        a5 = fitCoeffs(10);
        b5 = fitCoeffs(11);
        a6 = fitCoeffs(12);
        b6 = fitCoeffs(13);
        a7 = fitCoeffs(14);
        b7 = fitCoeffs(15);
        a8 = fitCoeffs(16);
        b8 = fitCoeffs(17);
        w = fitCoeffs(18);

        % builds up fit function
        fitFunctionBLanon = @(x) a0 + a1*cos(x*w) + b1*sin(x*w) + ...
            a2*cos(2*x*w) + b2*sin(2*x*w) + a3*cos(3*x*w) + b3*sin(3*x*w) + ...
            a4*cos(4*x*w) + b4*sin(4*x*w) + a5*cos(5*x*w) + b5*sin(5*x*w) + ...
            a6*cos(6*x*w) + b6*sin(6*x*w) + a7*cos(7*x*w) + b7*sin(7*x*w) + ...
            a8*cos(8*x*w) + b8*sin(8*x*w);
        fitFunction = @(x) fitFunctionBLanon(x)*0.5.*(tanh(-1e100*(x-n(end)))+1)+...
            quantityEnd*0.5.*(tanh(1e100*(x-n(end)))+1);

        aStringCoeffs = fitCoeffs(2:2:16);
        bStringCoeffs = fitCoeffs(3:2:17);

        % constructs string suitable for PyFR
        fitString = ['(',num2str(a0),'+'];
        for i = 1:7
            fitString = [fitString,num2str(aStringCoeffs(i)),'*cos(',num2str(w),'*',num2str(i),'*y)+',...
                num2str(bStringCoeffs(i)),'*sin(',num2str(w),'*',num2str(i),'*y)+'];
        end
        fitString = [fitString,num2str(aStringCoeffs(8)),'*cos(',num2str(w),'*',num2str(8),'*y)+',...
            num2str(bStringCoeffs(8)),'*sin(',num2str(w),'*',num2str(8),'*y)'];
        fitString = [fitString,')*0.5*(tanh(-10000000*(y-',num2str(n(end)),'))+1)+',...
            num2str(quantityEnd),'*0.5*(tanh(10000000*(y-',num2str(n(end)),'))+1)'];
    
    case 'sin8'
        fitFunctionBL = fit(n,quantityArray,'sin8','Normalize','on');

        % obtains coefficient values
        fitCoeffs = coeffvalues(fitFunctionBL);
        a1 = fitCoeffs(1);
        b1 = fitCoeffs(2);
        c1 = fitCoeffs(3);
        a2 = fitCoeffs(4);
        b2 = fitCoeffs(5);
        c2 = fitCoeffs(6);
        a3 = fitCoeffs(7);
        b3 = fitCoeffs(8);
        c3 = fitCoeffs(9);
        a4 = fitCoeffs(10);
        b4 = fitCoeffs(11);
        c4 = fitCoeffs(12);
        a5 = fitCoeffs(13);
        b5 = fitCoeffs(14);
        c5 = fitCoeffs(15);
        a6 = fitCoeffs(16);
        b6 = fitCoeffs(17);
        c6 = fitCoeffs(18);
        a7 = fitCoeffs(19);
        b7 = fitCoeffs(20);
        c7 = fitCoeffs(21);
        a8 = fitCoeffs(22);
        b8 = fitCoeffs(23);
        c8 = fitCoeffs(24);
        
        nMean = mean(n);
        nStd = std(n);

        % builds up fit function
        fitFunctionBLanon = @(x) a1*sin(b1*((x-nMean)/nStd)+c1) + ...
            a2*sin(b2*((x-nMean)/nStd)+c2) + a3*sin(b3*((x-nMean)/nStd)+c3) + ...
            a4*sin(b4*((x-nMean)/nStd)+c4) + a5*sin(b5*((x-nMean)/nStd)+c5) + ...
            a6*sin(b6*((x-nMean)/nStd)+c6) + ...
            a7*sin(b7*((x-nMean)/nStd)+c7) + a8*sin(b8*((x-nMean)/nStd)+c8);
        fitFunction = @(x) fitFunctionBLanon(x)*0.5.*(tanh(-1e100*(x-n(end)))+1)+...
            quantityEnd*0.5.*(tanh(1e100*(x-n(end)))+1);

        aStringCoeffs = fitCoeffs(1:3:22);
        bStringCoeffs = fitCoeffs(2:3:23);
        cStringCoeffs = fitCoeffs(3:3:24);

        % constructs string suitable for PyFR
        fitString = '(';
        for i = 1:7
            fitString = [fitString,num2str(aStringCoeffs(i)),'*sin(',num2str(bStringCoeffs(i)),...
                '*((y-',num2str(nMean),')/',num2str(nStd),')+',...
                num2str(cStringCoeffs(i)),')+'];
        end
        fitString = [fitString,num2str(aStringCoeffs(8)),'*sin(',num2str(bStringCoeffs(8)),...
                '*((y-',num2str(nMean),')/',num2str(nStd),')+',...
                num2str(cStringCoeffs(8)),')'];
        fitString = [fitString,')*0.5*(tanh(-10000000*(y-',num2str(n(end)),'))+1)+',...
            num2str(quantityEnd),'*0.5*(tanh(10000000*(y-',num2str(n(end)),'))+1)'];
    
    case 'generalLogistic'
        lowerBound = quantityArray(1);
        upperBound = quantityEnd;
        generalLogisticFitType = fittype( @(b,d,q,v,x) ...
            lowerBound+((upperBound-lowerBound)./((1+q*exp(-b*x+d)).^(1/v))));
        fitFunctionBL = fit(n,quantityArray,generalLogisticFitType,'Lower',[-inf -inf 0 -inf],'Upper',[inf inf inf inf],...
            'StartPoint',startPoints,'Display','iter','MaxFunEvals',5000,'MaxIter',5000);
        % obtains coefficient values
        fitCoeffs = coeffvalues(fitFunctionBL);
        bCoeff = fitCoeffs(1);
        dCoeff = fitCoeffs(2);
        qCoeff = fitCoeffs(3);
        vCoeff = fitCoeffs(4);

        fitFunctionBLanon = @(x) lowerBound+((upperBound-lowerBound)./((1+qCoeff*exp(-bCoeff*x+dCoeff)).^(1/vCoeff)));
        fitFunction = @(x) fitFunctionBLanon(x)*0.5.*(tanh(-1e100*(x-n(end)))+1)+...
            quantityEnd*0.5.*(tanh(1e100*(x-n(end)))+1);
        fitString = ['(',num2str(lowerBound),'+((',num2str(upperBound),'-',num2str(lowerBound),')/(pow((1+',num2str(qCoeff),...
            '*exp(-',num2str(bCoeff),'*y+',num2str(dCoeff),')),(1/',num2str(vCoeff),'))))',')*0.5*(tanh(-10000000*(y-',num2str(n(end)),'))+1)+',...
            num2str(quantityEnd),'*0.5*(tanh(10000000*(y-',num2str(n(end)),'))+1)'];

    case 'generalLogisticUpper'
        upperBound = quantityEnd;
        generalLogisticFitType = fittype( @(b,d,q,v,k,x) ...
            k+((upperBound-k)./((1+q*exp(-b*x+d)).^(1/v))));
        fitFunctionBL = fit(n,quantityArray,generalLogisticFitType,'Lower',[-inf -inf 0 -inf -inf],'Upper',[inf inf inf inf inf],...
            'StartPoint',startPoints,'Display','iter','MaxFunEvals',5000,'MaxIter',5000);
        % obtains coefficient values
        fitCoeffs = coeffvalues(fitFunctionBL);
        bCoeff = fitCoeffs(1);
        dCoeff = fitCoeffs(2);
        qCoeff = fitCoeffs(3);
        vCoeff = fitCoeffs(4);
        kCoeff = fitCoeffs(5);

        fitFunctionBLanon = @(x) kCoeff+((upperBound-kCoeff)./((1+qCoeff*exp(-bCoeff*x+dCoeff)).^(1/vCoeff)));
        fitFunction = @(x) fitFunctionBLanon(x)*0.5.*(tanh(-1e100*(x-n(end)))+1)+...
            quantityEnd*0.5.*(tanh(1e100*(x-n(end)))+1);
        fitString = ['(',num2str(kCoeff),'+((',num2str(upperBound),'-(',num2str(kCoeff),'))/(pow((1+',num2str(qCoeff),...
            '*exp(-',num2str(bCoeff),'*y+',num2str(dCoeff),')),(1/',num2str(vCoeff),'))))',')*0.5*(tanh(-10000000*(y-',num2str(n(end)),'))+1)+',...
            num2str(quantityEnd),'*0.5*(tanh(10000000*(y-',num2str(n(end)),'))+1)'];
end
