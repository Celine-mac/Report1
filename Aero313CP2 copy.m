clc;
close all;

nu=1.0;
h=1.0;
Ti=10.0;
Q1=10.0;
Q2_ref=5.0;
n=12;

A=[ 2 -1 0 0;
    1 2 -1 0;
    0 1 2 -1;
    0 0 1 2];

d = [0.0 6.0 2.0 -5.0];
[a,b,c] = LUdecomp(A);
[x] = LUsolve(a,b,c,d);

dts= [0.01, 0.10, 0.50];
[y1, t1, u1]= heatSolver(0.1, 0.01, Q2_ref);
[y2, t2, u2]= heatSolver(0.1, 0.10, Q2_ref);
[y3, t3, u3]= heatSolver(0.1, 0.50, Q2_ref);

Ns= [5,10,20,30,60,120,200];
for i=1:length(Ns)
    t=Ns(i)*0.01;
    plot(y1, anaSol(y1, t , n, Q2_ref),'--')
    hold on 
%    plot(y1, u1(Ns(i),:), 'k-')
end
xlabel("y");
ylabel("Temperature")

figure()
subplot(2,1,1)
plot(t1, u1(:,1),'-')
hold on
subplot(2,3,5)
plot(t2,u2(:,1),'-')
plot(t3,u3(:,1))
plot(t1,anaSol(0.0,t1,n,Q2_ref),'-')
ylabel("Temperature")

subplot(2,1,2)
plot(t1, u1(:,1)- anaSol(0.0,t1,n,Q2_ref),'-')
hold on
plot(t2,u2(:,1)-anaSol(0.0,t2,n,Q2_ref),'-')
plot(t3,u3(:,1)-anaSol(0.0,t3,n,Q2_ref),'-')

Qs=linspace(5.0, 15.0, 31);
us=[];

for i=1:length(Qs)
    [~, ~, u]=heatSolver(0.1, 0.5, Qs(i));
    us= [us, u(length(u),1)];
end
figure()
plot(Qs, us, 'b-', "DisplayName", "Numericall Solution")
xlabel("Q")
ylabel("max temperture")

%Function definitions
function [a, b, c] = LUdecomp(T)
 
a=diag(T);
b=diag(T,-1);
c=diag(T,1);
N=length(a);

for i=2:N
    b(i-1)=b(i-1)/a(i-1);
    a(i)=a(i)-b(i-1)*c(i-1);
end

end

function[x] = LUsolve(a,b,c,d)
x=zeros(size(d),'like',d);
N=length(d);

x(1)=d(1);

for j=2:N
    x(j)=d(j)-b(j-1)*x(j-1);
    x(N)=x(N)/a(N);
end
for k = N-1:-1:1
    x(k)=(x(k)-c(k)*x(k+1))/a(k);
end

end

function[u]=anaSol(y,t,n,Q2)

nu=1.0;
h=1.0;
Q1=10.0;

u=zeros(size(y),'like',y);
w=(Q1-Q2)/(2*h)*y.^2-Q1*y;
u=u+w;
u=u+(2*Q1+Q2)*h/6;
for i=1:n+1
    p=i*pi/h;
    if mod(i,2)==0
        An = (2*h*(Q2-Q1))/(i*pi)^2;
    else
        An = (2*h*(-Q1-Q2))/ (i*pi)^2;
    end
    u=u+An*cos(p*y)*exp(-nu*p^2*t);
end
u=u+nu*(Q1-Q2)/h*t;
end

    function[y, t, U]=heatSolver(dy, dt, Q2)

h=1.0;
Ti=10.0;

Ny=int16(ceil(h/dy));
Nt=int16(ceil(Ti/dt));
y=linspace(0,h,Ny+1);
t=linspace(0,Ti,Nt+1);
U=zeros(Nt+1, Ny+1);

[F, T, S, B1, B2]= initProblem(Ny, dy, dt, Q2);
[a,b,c]= LUdecomp(T);

for i=1:Nt
    u=U(i, 2:Ny)';
    U(i+1, 2:Ny) = LUsolve(a,b,c, S*u+F);
    U(i+1, 1)=B1(1)+ B1(2)*U(i+1,2) + B1(3)*U(i+1, 3);
    U(i+1, Ny+1)= B2(1)+B2(2)*U(i+1, length(U(i+1, :))-1) +B2(3)*U(i+1, length(U(i+1, :))-2);



end


    end

    function [F, T, S, B1, B2]= initProblem(Ny, dy, dt, Q2)
        Q1= 10.0;
        nu= 1.0;

        N= Ny -1;
        f=zeros(N, 1);
        f(1)= Q1*2.0/3.0;
        f(N)= -Q2*2.0/3.0;

        A= -2.0*diag(ones(1,nu),0)+ diag(ones(1,N-1),1)+diag(ones(1, N-1),-1);
        A0= zeros(1,N);
        A0(2)= -A0(1);
        A(1,:)= A0;
        A(N,:)=fliplr(A0);
        I=diag(ones(1,N),0);
        r=(nu/2)*(dt/dy^2);

        F=dt/dy*nu*f;
        T=I-r*A;
        S=I+r*A;

        B1=[Q1*(2*dy)/3.0 4.0/3.0 -1.0/3.0];
        B2=[-Q1*(2*dy)/3.0 4.0/3.0 -1.0/3.0];

    end


