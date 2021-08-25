load('Decessi_ISTAT_comuni.mat');


%https://it.mathworks.com/help/matlab/ref/fprintf.html
%https://it.mathworks.com/help/matlab/ref/fopen.html


%%
%dataset

file_name='info.txt';
file_id=fopen(file_name,'w');

t1="pro_com";
t2="comune";
t3="cod_reg";
t4="cod_prov";
t5="tot_pop";
t6="density";
t7="x";
t8="y";

title=[t1 t2 t3 t4 t5 t6 t7 t8];

formatspec_title='%s,%s,%s,%s,%s,%s,%s,%s\r\n';

fprintf(file_id,formatspec_title,title);

format_spec='%u,%s,%u,%u,%u,%f,%e,%e\r\n';

for iii=1:7903
    pro_com=Record(iii).PRO_COM;
    comune=Record(iii).COMUNE;
    cod_reg=Record(iii).COD_REG;
    cod_prov=Record(iii).COD_PROV;
    tot_pop=Record(iii).Tot_Pop;
    density=Record(iii).Densita;
    x=Record(iii).Xc;
    y=Record(iii).Yc;
    
    line=[cod_reg cod_prov tot_pop density x y];
    fprintf(file_id,format_spec,pro_com,comune,line);
    
end

fclose(file_id);


%% indexes


file_name='indexes.txt';
file_id=fopen(file_name,'w');
 fprintf(file_id,'%s\n','indexes');
for i=1:length(indexes_local)
    fprintf(file_id,'%u\n',indexes_local(1,i));
end

fclose(file_id);

%% dataset-covariates
load('Decessi_ISTAT_comuni.mat');
load('coordinates_from_static_dataset.mat')

file_name='dataset_covariates.txt';
file_id=fopen(file_name,'w');

t1="id";%pro_com
t2="name";%comune
t3="cod_reg";

t4="x";
t5="y";

t6="population";
t7="density";

t8="mean_norm_deaths";
t9="mean_norm_deaths_covid";

title=[t1 t2 t3 t4 t5 t6 t7 t8,t9];

formatspec_title='%s,%s,%s,%s,%s,%s,%s,%s,%s\r\n';

fprintf(file_id,formatspec_title,title);

format_spec='%u,%s,%u,%.6f,%.6f,%u,%.6f,%.9f,%.9f\r\n';

for iii=1:7903
    pro_com=Record(iii).PRO_COM;
    comune=Record(iii).COMUNE;
    cod_reg=Record(iii).COD_REG;
    x=coordinatesfromstaticDS(iii,1);
    y=coordinatesfromstaticDS(iii,2);
    tot_pop=Record(iii).Tot_Pop;
    density=Record(iii).Densita;

    M=Record(iii).Tot_Time;
    mnd=sum(1/9* ones(1,9)*M(1:9,:))/tot_pop;
    mndc=sum(M(10,:))/tot_pop;

    fprintf(file_id,format_spec,pro_com,comune,cod_reg, x, y, tot_pop, density, mnd, mndc);
    
end

fclose(file_id);


