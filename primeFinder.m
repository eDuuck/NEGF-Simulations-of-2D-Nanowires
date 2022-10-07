function [primeVec,elapsedTime] = primeFinder(limit)
    tic
    if limit == 1
        primeVec = [];
        elapsedTime = toc;
        return
    end

    primeVec = zeros(1,ceil((limit/log(limit))*(1 + 1.2762/log(limit))));
    
    primeList = load("primes.mat").primes;
    primes = length(primeList);
    primeVec(1:primes) = primeList;

    if limit <= primeVec(primes)
        idx = find(primeVec <= limit,1,"last");
        primeVec = primeVec(1:idx);
        elapsedTime = toc;
        return
    end
    currLim = floor(sqrt(primeVec(primes)));
    limUpdate = (currLim+1)^2;
    newTestVec = false;

    idx = find(primeVec <= currLim,1,"first");
    testVec = primeVec(1:idx);
    
    testList = primeVec(primes)+1:limit;
    for i = 1:primes
        testList(mod(testList,primeVec(i)) == 0) = [];
    end
    for n = testList

        if limUpdate <= n
            currLim = currLim + 1;
            limUpdate = (currLim+1)^2;
            newTestVec = true;
        end
        testLim = currLim;

        if newTestVec
            if primeVec(idx+1) <= testLim
                idx = idx+1;
                testVec = primeVec(1:idx);
            end
        end
        isprime = true;
        for d = testVec
            if mod(n,d) == 0
                isprime = false;
                break;
            end
        end
        if isprime
            primes = primes+1;
            primeVec(primes) = n;
        end
    end
    elapsedTime = toc;
    primeVec = primeVec(1:primes);
end