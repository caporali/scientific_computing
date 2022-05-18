function plot_function(alpha,n,err,k,t,x)
% Support function to plot GR estimated solutions
%
% 	input:
%       alpha := damping factor
%       n := n. of the method
% 		err := error array 
% 		k := n. matrix-vector products
%       t := time
%       x := estimate solution 
% 
% 	example
%		plot_function(alpha,1,err,k,t,x)

if n == 1
    % plot
    semilogy(err,"--");

    title(['Errore al passo n-esimo / $$\alpha = ',num2str(alpha),'$$'],'interpreter','latex');
    xlabel('n. matrix-vector products','interpreter','latex');
    ylabel('$$||r||_1$$','interpreter','latex');

    ax = gca;
    ax.XAxis.FontSize = 8;
    ax.YAxis.FontSize = 8;
    ax.XLabel.FontSize = 10;
    ax.YLabel.FontSize = 10;
    ax.Title.FontSize = 12;

    hold on;

    % text output: time, number of matrix-vector products and x
    fprintf("\n	n. matrix-vector products = %d \n", k);
    fprintf("	t = %f s \n", toc);
    fprintf("	x(1:3,1) = [%e %e %e ...]' \n", x(1), x(2), x(3));

else
    % plot
    if (n == 2) || (n == 3) || (n == 4)        
        semilogy(err,"--");
    else
        semilogy(err,"-");
    end
    hold on;

    % text output: time, number of matrix-vector products and x
    fprintf("\n	n. matrix-vector products = %d \n", k);
    fprintf("	t = %f s \n", t);
    fprintf("	x(1:3,1) = [%e %e %e ...]' \n", x(1), x(2), x(3));
end