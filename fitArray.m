function [fitFunction, fitString] = fitArray(n,quantityArray,quantityEnd,fitType,startPoints,nSplit)
% quantity end is rhoE, vE, uE etc
% fitType can be fourier8, sin8, general logistic, general logistic only
% upper bound
% nSplit is only used for split fit types

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

    case 'splitPolySqrt'
        % first doing polynomial fit for lower section
        % splitting off lower section
        nLower = n(n<nSplit);
        quantityArrayLower = quantityArray(n<nSplit);
        fitFunctionBLlower = fit(nLower,quantityArrayLower,'poly8','Normalize','on');
        nLowerMean = mean(nLower);
        nLowerStd = std(nLower);
        
        % obtains coefficient values for lower
        fitCoeffsLower = coeffvalues(fitFunctionBLlower);
        p1Coeff = fitCoeffsLower(1);
        p2Coeff = fitCoeffsLower(2);
        p3Coeff = fitCoeffsLower(3);
        p4Coeff = fitCoeffsLower(4);
        p5Coeff = fitCoeffsLower(5);
        p6Coeff = fitCoeffsLower(6);
        p7Coeff = fitCoeffsLower(7);
        p8Coeff = fitCoeffsLower(8);
        p9Coeff = fitCoeffsLower(9);

        fitFunctionBLlowerAnon = @(x) p1Coeff*((x-nLowerMean)./nLowerStd).^8+...
        p2Coeff*((x-nLowerMean)./nLowerStd).^7+...
        p3Coeff*((x-nLowerMean)./nLowerStd).^6+...
        p4Coeff*((x-nLowerMean)./nLowerStd).^5+... 
        p5Coeff*((x-nLowerMean)./nLowerStd).^4+...
        p6Coeff*((x-nLowerMean)./nLowerStd).^3+...
        p7Coeff*((x-nLowerMean)./nLowerStd).^2+...
        p8Coeff*((x-nLowerMean)./nLowerStd)+...
        p9Coeff;

        % sqrt root fit for upper section
        nUpper = n(not(n<nSplit));
        quantityArrayUpper = quantityArray(not(n<nSplit));
        sqrtFitType = fittype( @(a,b,c,d,x) ...
            c*x-sqrt((c*x+b).^2+a)+d );
        fitFunctionBLupper = fit(nUpper,quantityArrayUpper,sqrtFitType,'Lower',[0 -inf -inf -inf],...
            'StartPoint',startPoints,'Display','iter','MaxFunEvals',5000,'MaxIter',5000);

        fitCoeffsUpper = coeffvalues(fitFunctionBLupper);
        aCoeff = fitCoeffsUpper(1);
        bCoeff = fitCoeffsUpper(2);
        cCoeff = fitCoeffsUpper(3);
        dCoeff = fitCoeffsUpper(4);

        fitFunctionBLupperAnon = @(x) cCoeff*x-sqrt((cCoeff*x+bCoeff).^2+aCoeff)+dCoeff;

        % assembling overall function
        fitFunction = @(x) fitFunctionBLlowerAnon(x).*0.5.*(tanh(-1e100*(x-nSplit))+1)+...
            fitFunctionBLupperAnon(x).*(0.5*(tanh(1e100*(x-nSplit))+1) + 0.5*(tanh(-1e100*(x-n(end)))+1) -1)+...
            quantityEnd*0.5.*(tanh(1e100*(x-n(end)))+1);

        % generating PyFR string
        fitStringLower = ['(',num2str(p1Coeff),'*pow(((y-',num2str(nLowerMean),')/',num2str(nLowerStd),'),8)+',...
            num2str(p2Coeff),'*pow(((y-',num2str(nLowerMean),')/',num2str(nLowerStd),'),7)+',...
            num2str(p3Coeff),'*pow(((y-',num2str(nLowerMean),')/',num2str(nLowerStd),'),6)+',...
            num2str(p4Coeff),'*pow(((y-',num2str(nLowerMean),')/',num2str(nLowerStd),'),5)+',...
            num2str(p5Coeff),'*pow(((y-',num2str(nLowerMean),')/',num2str(nLowerStd),'),4)+',...
            num2str(p6Coeff),'*pow(((y-',num2str(nLowerMean),')/',num2str(nLowerStd),'),3)+',...
            num2str(p7Coeff),'*pow(((y-',num2str(nLowerMean),')/',num2str(nLowerStd),'),2)+',...
            num2str(p8Coeff),'*((y-',num2str(nLowerMean),')/',num2str(nLowerStd),')+',num2str(p9Coeff),...
            ')*0.5*(tanh(-10000000*(y-',num2str(nSplit),'))+1)'];
        fitStringUpper = ['(',num2str(cCoeff),'*y-sqrt(pow((',num2str(cCoeff),'*y+',num2str(bCoeff),...
            '),2)+',num2str(aCoeff),')+',num2str(dCoeff),')*(0.5*(tanh(10000000*(y-',num2str(nSplit),...
            '))+1)+0.5*(tanh(-10000000*(y-',num2str(n(end)),'))+1)-1)'];
        fitString = [fitStringLower,'+',fitStringUpper,'+',...
            num2str(quantityEnd),'*0.5*(tanh(10000000*(y-',num2str(n(end)),'))+1)'];
    
    case 'splitPolyGeneralLogistic'

        % first doing polynomial fit for lower section
        % splitting off lower section
        nLower = n(n<nSplit);
        quantityArrayLower = quantityArray(n<nSplit);
        fitFunctionBLlower = fit(nLower,quantityArrayLower,'poly8','Normalize','on');
        nLowerMean = mean(nLower);
        nLowerStd = std(nLower);
        
        % obtains coefficient values for lower
        fitCoeffsLower = coeffvalues(fitFunctionBLlower);
        p1Coeff = fitCoeffsLower(1);
        p2Coeff = fitCoeffsLower(2);
        p3Coeff = fitCoeffsLower(3);
        p4Coeff = fitCoeffsLower(4);
        p5Coeff = fitCoeffsLower(5);
        p6Coeff = fitCoeffsLower(6);
        p7Coeff = fitCoeffsLower(7);
        p8Coeff = fitCoeffsLower(8);
        p9Coeff = fitCoeffsLower(9);

        fitFunctionBLlowerAnon = @(x) p1Coeff*((x-nLowerMean)./nLowerStd).^8+...
        p2Coeff*((x-nLowerMean)./nLowerStd).^7+...
        p3Coeff*((x-nLowerMean)./nLowerStd).^6+...
        p4Coeff*((x-nLowerMean)./nLowerStd).^5+... 
        p5Coeff*((x-nLowerMean)./nLowerStd).^4+...
        p6Coeff*((x-nLowerMean)./nLowerStd).^3+...
        p7Coeff*((x-nLowerMean)./nLowerStd).^2+...
        p8Coeff*((x-nLowerMean)./nLowerStd)+...
        p9Coeff;

        % general logistic for upper section
        nUpper = n(not(n<nSplit));
        quantityArrayUpper = quantityArray(not(n<nSplit));

        upperBound = quantityEnd;
        generalLogisticFitType = fittype( @(b,d,q,v,k,x) ...
            k+((upperBound-k)./((1+q*exp(-b*x+d)).^(1/v))));
        fitFunctionBLupper = fit(nUpper,quantityArrayUpper,generalLogisticFitType,'Lower',[-inf -inf -inf -inf 0],'Upper',[inf inf inf inf inf],...
            'StartPoint',startPoints,'Display','iter','MaxFunEvals',5000,'MaxIter',5000);
        % obtains coefficient values
        fitCoeffs = coeffvalues(fitFunctionBLupper);
        bCoeff = fitCoeffs(1);
        dCoeff = fitCoeffs(2);
        qCoeff = fitCoeffs(3);
        vCoeff = fitCoeffs(4);
        kCoeff = fitCoeffs(5);

        fitFunctionBLupperAnon = @(x) kCoeff+((upperBound-kCoeff)./((1+qCoeff*exp(-bCoeff*x+dCoeff)).^(1/vCoeff)));

        % assembling overall function
        fitFunction = @(x) fitFunctionBLlowerAnon(x).*0.5.*(tanh(-1e100*(x-nSplit))+1)+...
            fitFunctionBLupperAnon(x).*(0.5*(tanh(1e100*(x-nSplit))+1) + 0.5*(tanh(-1e100*(x-n(end)))+1) -1)+...
            quantityEnd*0.5.*(tanh(1e100*(x-n(end)))+1);

        % generating PyFR string
        fitStringLower = ['(',num2str(p1Coeff),'*pow(((y-',num2str(nLowerMean),')/',num2str(nLowerStd),'),8)+',...
            num2str(p2Coeff),'*pow(((y-',num2str(nLowerMean),')/',num2str(nLowerStd),'),7)+',...
            num2str(p3Coeff),'*pow(((y-',num2str(nLowerMean),')/',num2str(nLowerStd),'),6)+',...
            num2str(p4Coeff),'*pow(((y-',num2str(nLowerMean),')/',num2str(nLowerStd),'),5)+',...
            num2str(p5Coeff),'*pow(((y-',num2str(nLowerMean),')/',num2str(nLowerStd),'),4)+',...
            num2str(p6Coeff),'*pow(((y-',num2str(nLowerMean),')/',num2str(nLowerStd),'),3)+',...
            num2str(p7Coeff),'*pow(((y-',num2str(nLowerMean),')/',num2str(nLowerStd),'),2)+',...
            num2str(p8Coeff),'*((y-',num2str(nLowerMean),')/',num2str(nLowerStd),')+',num2str(p9Coeff),...
            ')*0.5*(tanh(-10000000*(y-',num2str(nSplit),'))+1)'];

        fitStringUpper = ['(',num2str(kCoeff),'+((',num2str(upperBound),'-(',num2str(kCoeff),'))/(pow((1+',num2str(qCoeff),...
            '*exp(-',num2str(bCoeff),'*y+',num2str(dCoeff),')),(1/',num2str(vCoeff),...
            '))))',')*(0.5*(tanh(10000000*(y-',num2str(nSplit),...
            '))+1)+0.5*(tanh(-10000000*(y-',num2str(n(end)),'))+1)-1)'];
        
        fitString = [fitStringLower,'+',fitStringUpper,'+',...
            num2str(quantityEnd),'*0.5*(tanh(10000000*(y-',num2str(n(end)),'))+1)'];
end
