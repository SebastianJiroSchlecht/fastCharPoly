%% La Budde's Method For Computing Characteristic Polynomials
% by Rizwana Rehman and Ilse C.F. Ipsen
% in arXiv, April 2011
%
% *Reproduce in Code*
% (c) Sebastian Jiro Schlecht:  Monday, 10. December 2018
%
clear; clc; close all;

%% Introduction
% In this example, we reproduce a fast algorithm (Algorithm 2) to computer the
% characteristic polynomial of a square matrix and compare the output to the 
% MATLAB build-in function charPoly.

N = 12;
A = randn(N) + 1i*randn(N);
A = A / N;
disp('Comparison of Matlab charpoly and LaBudde charpoly')

%%-----------------
disp('Matlab')
tic
pMatlab = charpoly(A);
toc
disp(pMatlab)

%%-----------------
disp('LaBudde')
tic
pFast = fastCharPoly( A );
toc
disp(pFast)

%%compare
disp('Compare error')
disp(pMatlab - pFast)


%% Speed comparison
% We compare the performance of fastCharPoly to the MATLAB build-in function charPoly.
% It can be seen that fastCharPoly is more than 100 times faster.

matrixSizes = 1:16;
numberOfRepetitions = 5;
performanceTime_MATLAB = zeros(length(matrixSizes),numberOfRepetitions);
performanceTime_FAST = zeros(length(matrixSizes),numberOfRepetitions);
for N = matrixSizes
    % disp(matrixSize)
    for repetition = 1:numberOfRepetitions
        A = randn(N) + 1i*randn(N);
        A = A / N;
        tic
            pMatlab = charpoly( A );
        performanceTime_MATLAB(N,repetition) = toc;
        tic
            pFast = fastCharPoly( A );
        performanceTime_FAST(N,repetition) = toc;
        
        if max(abs(pMatlab - pFast)./pMatlab) < 1e-7
           % disp('OK')
        else
            warning(['Results deviate by ' num2str(max(abs(pMatlab - pFast)))])
        end
    end
end
    

%% Plot performance results
averagePerformanceTime_MATLAB = mean(performanceTime_MATLAB,2);
averagePerformanceTime_FAST = mean(performanceTime_FAST,2);

figure(1); hold on; grid on;
plot(matrixSizes,averagePerformanceTime_MATLAB)
plot(matrixSizes,averagePerformanceTime_FAST)
legend({'MATLAB','FAST'})
xlabel('Matrix Size')
ylabel('Time [s]')
set(gca,'YScale','log');
axis tight;