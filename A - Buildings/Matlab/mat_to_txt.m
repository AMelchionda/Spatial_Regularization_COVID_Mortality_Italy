load('Decessi_ISTAT_comuni.mat')

main = cell2table(num2cell(Record(1).Tot_Time(:,:)));

C = {Record(1).PRO_COM, Record(1).PRO_COM, Record(1).PRO_COM, ...
    Record(1).PRO_COM, Record(1).PRO_COM, Record(1).PRO_COM, ...
    Record(1).PRO_COM, Record(1).PRO_COM, Record(1).PRO_COM, Record(1).PRO_COM};

D = 2011:2020;

E = {Record(1).Tot_Pop, Record(1).Tot_Pop, Record(1).Tot_Pop, ...
    Record(1).Tot_Pop, Record(1).Tot_Pop, Record(1).Tot_Pop, ...
    Record(1).Tot_Pop, Record(1).Tot_Pop, Record(1).Tot_Pop, Record(1).Tot_Pop};

F = {Record(1).Xc, Record(1).Xc, Record(1).Xc, ...
    Record(1).Xc, Record(1).Xc, Record(1).Xc, ...
    Record(1).Xc, Record(1).Xc, Record(1).Xc, Record(1).Xc};

G = {Record(1).Yc, Record(1).Yc, Record(1).Yc, ...
    Record(1).Yc, Record(1).Yc, Record(1).Yc, ...
    Record(1).Yc, Record(1).Yc, Record(1).Yc, Record(1).Yc};

main.COMUNE = C';

main.YEAR = D';

main.TOT_POP = E';

main.COORD_X = F';

main.COORD_Y = G';

main = [main(:,end) main(:,end-1) main(:,end-2) main(:,end-3) main(:,end-4) main(:,1:end-5)];

i = 2;

while i < 7904

    T = cell2table(num2cell(Record(i).Tot_Time(:,:)));

    C = {Record(i).PRO_COM, Record(i).PRO_COM, Record(i).PRO_COM, ...
    Record(i).PRO_COM, Record(i).PRO_COM, Record(i).PRO_COM, ...
    Record(i).PRO_COM, Record(i).PRO_COM, Record(i).PRO_COM, Record(i).PRO_COM};

    E = {Record(i).Tot_Pop, Record(i).Tot_Pop, Record(i).Tot_Pop, ...
    Record(i).Tot_Pop, Record(i).Tot_Pop, Record(i).Tot_Pop, ...
    Record(i).Tot_Pop, Record(i).Tot_Pop, Record(i).Tot_Pop, Record(i).Tot_Pop};
    
    F = {Record(i).Xc, Record(i).Xc, Record(i).Xc, ...
    Record(i).Xc, Record(i).Xc, Record(i).Xc, ...
    Record(i).Xc, Record(i).Xc, Record(i).Xc, Record(i).Xc};

    G = {Record(i).Yc, Record(i).Yc, Record(i).Yc, ...
    Record(i).Yc, Record(i).Yc, Record(i).Yc, ...
    Record(i).Yc, Record(i).Yc, Record(i).Yc, Record(i).Yc};

    T.COMUNE = C';

    T.YEAR = D';

    T.TOT_POP = E';

    T.COORD_X = F';

    T.COORD_Y = G'; 

    T = [T(:,end) T(:,end-1) T(:,end-2) T(:,end-3) T(:,end-4) T(:,1:end-5)];

    main = [main; T];

    i = i+1;

    disp(i)

end

writetable(main,'Complete.txt')