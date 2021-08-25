
%%

%EPSG23032

load('Italy_boundary.mat')
load('sardegna_boundary')
load('coordinates_rescaled_mat.mat')


X=[coordinatesrescaled(:,1);italyboundary(:,1);sardegnaboundary(:,1)];
Y=[coordinatesrescaled(:,2);italyboundary(:,2);sardegnaboundary(:,2)];

N=length(X);

figure()
disp('EPSG23032 (originale)')
%%

%EPSG23032

load('Italy_boundary.mat')
load('sardegna_boundary')
load('coordinates_rescaled_mat.mat')


X=[coordinatesrescaled(:,1)];
Y=[coordinatesrescaled(:,2)];

N=length(X);

figure()
disp('EPSG23032 (originale)')


%%
load('mainland_boundary_new.mat')
load('island_boundary_new')
load('coordinates_rescaled_mat.mat')


X=[coordinatesrescaled(:,1)];
Y=[coordinatesrescaled(:,2)];

italyboundary=mainlandboundarynew;
sardegnaboundary=islandboundarynew;

N=length(X);

figure()
disp('EPSG23032 (originale)')

%%
flag=true;
flag2=false;
indexes_global=[];

plot(X,Y,'k*')
hold on

plot(italyboundary(:,1),italyboundary(:,2),'r-')
plot(sardegnaboundary(:,1),sardegnaboundary(:,2),'r-')

indexes_local=[];

start=menu('Start?','yes','no');
if(start==1)
  while(flag)

        [x,y]=ginput(1);

        index=1;
        for i=2:N
            if (((X(i)-x)^2+(Y(i)-y)^2)<((X(index)-x)^2+(Y(index)-y)^2))
                index=i;
            end
        end
        
        
        if(sum(indexes_local==index)>0)
            cont=menu('Continue?','yes','no','delete points');
            if(cont==2)
                flag=false;
            else
                if(cont==3)
                    flag2=true;
                end
            end
            while(flag2)
                [x,y]=ginput(1);
                
                index=1;
                for i=2:N
                    if (((X(i)-x)^2+(Y(i)-y)^2)<((X(index)-x)^2+(Y(index)-y)^2))
                        index=i;
                    end
                end
                
                plot(X(index),Y(index),'xr','LineWidth',1.5)
                
                if(sum(indexes_local==index)>0)
                    v=1:length(indexes_local);
                    pos=v(indexes_local==index);
                    if(pos>1 && pos<v(end))
                        indexes_local=[indexes_local(1:pos-1) indexes_local(pos+1:end)];
                    end
                    if(pos==1)
                        indexes_local=indexes_local(2:end);
                    else
                        if(pos==v(end))
                            indexes_local=indexes_local(1:end-1);
                        end
                    end
                end
                cont2=menu('Continue deleting?','yes','no');
                if(cont2==2)
                    flag2=false;
                end
            end
        else
            indexes_local=[indexes_local,index];
            plot(X(indexes_local),Y(indexes_local),'-og','LineWidth',1.5)
        end

    end
end

boundary_pt=[X(indexes_local),Y(indexes_local)];

duplicated_nodes=false;
for i=1:length(indexes_global)
    if(sum(indexes_global==indexes_global(i))>1)
        duplicated_nodes=true;
    end
end

if(duplicated_nodes)
    disp("duplicated_nodes")
end

