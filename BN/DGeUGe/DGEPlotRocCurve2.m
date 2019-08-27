function DGEPlotRocCurve2(p_falseTotal,p_correctTotal,LimitX,n)

m=0:0.01:1;
p_falseTotal=fliplr(p_falseTotal);
p_correctTotal=fliplr(p_correctTotal);
if isempty(find(p_falseTotal==LimitX))
    %% is empty, no x=e then have to interpolate to insert a poinx(e,y)
    %% first put together x and y and tranpose it to sort
    c=cat(1,p_falseTotal,p_correctTotal)';
    % now sort in ascending order using x
    c_ord=sortrows(c,1);
    %find the value immediately after LimitX
    [a b]=find(c_ord>LimitX);
    %sum the y value immediately after and before LimitX
    %yLimitX=(c_ord(a(1),2)+c_ord(a(1)-1,2))/2;
    %Having two points that compose the line
    x1=c_ord(a(1)-1,1);
    x2=c_ord(a(1),1);
    y1=c_ord(a(1)-1,2);
    y2=c_ord(a(1),2);
    %and one point in x, LimitX, we can interpolate y
    yLimitX=(((LimitX-x1)/(x2-x1))*(y2-y1))+y1;
    %add a new point to the data
    p_falseTotal(end+1)=LimitX;
    p_correctTotal(end+1)=yLimitX;
    %p_falseTotal
    %p_correctTotal
    %% just get the points before e
    p_false=p_falseTotal(find(p_falseTotal<=LimitX));
    p_correct=p_correctTotal(find(p_falseTotal<=LimitX));
   
     
else
    %% just get the points before e
    p_false=p_falseTotal(find(p_falseTotal<=LimitX));
    p_correct=p_correctTotal(find(p_falseTotal<=LimitX));
   
end





figure(n)
plot(p_false,p_correct,'b-','Linewidth',2)

axis([-0.05*LimitX LimitX -0.05  1])
set (gca, 'Fontsize', 15)
ylabel('% True')
xlabel('% False')
area=abs(trapz(p_false,p_correct));
legend(['Area= ',num2str(area,5)])
hold on
plot(m,m,'b:','Linewidth',2)