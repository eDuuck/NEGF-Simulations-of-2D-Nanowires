method = 2;

points = 9;
tou = 0.05;
eta = 0.001;
E = 1.5;
Ec = 1;
eps = Ec + 4*tou;
I = eye(points);
gn = zeros(points);

alpha = eye(points)*eps - tou*(diag(ones(points-1,1),1) + diag(ones(points-1,1),-1));
beta = -eye(points)*tou;

iterations = 5000;
% g_inf = 0;
% 
% 
% 
% for j = 1:10000
%     g_inf = (E*I+i*eta*I-alpha-beta'*g_inf*beta)^-1;
% end
% g_inf
%%
eta = 0.001
rate = 0.7
converged = zeros(length(rate),length(eta));
time = zeros(length(rate),length(eta));
g_inf = zeros(points,points,length(eta));
m1 = zeros(1,length(eta));
m1t = zeros(1,length(eta));
acc = 1e-10;

for l = 1:length(eta)
    for j = 1:10000
        g_inf_new = (E*I+i*eta(l)*I-alpha-beta'*g_inf(:,:,l)*beta)^-1;
        if max(abs(g_inf_new - g_inf(:,:,l)),[],'all') < acc
            break
        end
        g_inf(:,:,l) = g_inf_new;
    end

    disp(l+"/"+length(eta));
    for k = 1:length(rate)
        %disp("Hey");
        %max_change = zeros(iterations,1);
        offset = zeros(iterations,1);
        gn = 0;
        tic
        for j = 1:iterations
            gnew = (E*I+i*eta(l)*I-alpha-beta'*gn*beta)^-1;
            change = gnew - gn;
            max_change(j) = max(abs(change),[],'all');
            if method == 1
                gn = gnew;
            elseif method == 2
                gn = gn + rate(k)*change;
            end
            offset(j) = sum(abs((g_inf(:,:,l) - gn).^2),'all')/points^2;
            if converged(k,l) == 0
                if max_change(j) < acc
                    converged(k,l) = j;
                    time(k,l) = toc;
                    break;
                end
            end
        end
        if converged(k,l) == 0
            converged(k,l) = 0;
            time(k,l) = toc;
        end
    end
end
%%
figure(1)
imagesc(abs(gn))
%%
% time = (converged~=0).*time;
% time = time + (converged==0).*(max(time,[],'all'));
% figure(2)
% subplot(2,2,1)
% imagesc('XData',eta,'YData',rate,'CData',log(time))
% title("Log of calculation time")
% subplot(2,2,2)
% imagesc('XData',eta,'YData',rate,'CData',log(converged))
% title("Log of converge iterations")
% 
% temp = ones(length(rate),1) * time(end,:) - time;
% temp = temp ./ (ones(length(rate),1) * time(end,:));
% temp = (temp > -0.1) .* temp;
% subplot(2,2,3)
% imagesc('XData',eta,'YData',rate,'CData',(temp))
% title("Ratio of saved time compared to simple algoritm");
% 
% temp = ones(length(rate),1) * converged(end,:) - converged;
% temp = (temp > 0) .* temp;
% subplot(2,2,4)
% imagesc('XData',eta,'YData',rate,'CData',log(temp))
% title("Log of saved iterations");
% 
% 
% for j = 1:4
%     subplot(2,2,j)
%     axis([eta(1) eta(end) rate(1) rate(end)]);
%     xlabel("Value of \eta");
%     ylabel("Rate of iteration change");
% end
