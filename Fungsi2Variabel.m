function GA_2

disp ('Algoritma Genetika : Fungsi dengan 2 variable')
disp (' ')
disp ('Masalah : Mencari nilai maksimum dari fungsi dengan 2 variable (1-x).^2.*exp(-x.^2-(y+1).^2-(x-x.^3.^5).*exp(-x.^2-y.^2)')
disp ('dimana parameter "x" dan "y" berada pada range -3 dan +3')

% Inisialisasi Populasi
jumlah_individu = 6;
jumlah_variabel = 2;
jumlah_gen = 16;
xymin = -3;
xymax = 3;
jumlah_generasi = 100;
probabilitas_crossover = 0.9;
probabilitas_mutasi = 0.005;


% pupulasi awal
chromosome = round(rand(jumlah_individu, jumlah_gen));

% menghitung nilai x (konversi biner ke desimal)
xy=zeros(jumlah_individu,jumlah_variabel);
jumlah_bit_variabel=jumlah_gen/jumlah_variabel;
if length(xymin)==1,xymin=xymin*ones(1,jumlah_variabel);end
if length(xymax)==1,xymax=xymax*ones(1,jumlah_variabel);end

for individu=1:jumlah_variabel
	xy(:,individu)=chromosome(:,1+jumlah_bit_variabel*(individu-1):jumlah_bit_variabel*individu)*[2.^(jumlah_bit_variabel-1:-1:0)]';
	xy(:,individu)=xymin(individu)+(xymax(individu)-xymin(individu))*(xy(:,individu)+1)./(2^jumlah_bit_variabel+1);
end

fungsi='(1-x).^2.*exp(-x.^2-(y+1).^2)-(x-x.^3-y.^5).*exp(-x.^2-y.^2)';

% Menghitung nilai fx
fx=evalFx(fungsi,xy(:,1),xy(:,2));
best(1)=max(fx);
ave(1)=mean(fx);

figure('name',['Lokasi Chromosome']);
[x,y]=meshgrid(xymin(1):.25:xymax(1),xymin(2):.25:xymax(2));
z=evalFx(fungsi,x,y); z=z+4;
mesh(x,y,z)
axis([-3 3 -3 3 0 6])
hold;
contour(x,y,z,20,'k')

scatter3(xy(:,1),xy(:,2),fx+4.08,40,'r','filled')
plot(xy(:,1),xy(:,2),'k.','markersize',23)
title(['Hit any key to continue']);
xlabel('Parameter "x"');
ylabel('Parameter "y"');
zlabel('Nilai Fitness');
hold;

for i = 1 : jumlah_generasi

	% Nilai Fitness
	nilai_fitness = fx;
	if min(fx) < 0
		nilai_fitness = fx - min(fx);
    end
    
% Seleksi
jumlah_individu_terbaik = round(jumlah_individu * probabilitas_crossover);
cumulative_fitness = repmat(cumsum(nilai_fitness), 1, jumlah_individu_terbaik);
chance = repmat(rand(1, jumlah_individu_terbaik), jumlah_individu, 1) * cumulative_fitness(jumlah_individu, 1);
[selind,j] = find(chance < cumulative_fitness & chance >= [zeros(1, jumlah_individu_terbaik); cumulative_fitness(1: jumlah_individu-1,:)]);
chromosome_baru = chromosome(selind,:);

% Crossover
point_crossover = round(rand(floor(jumlah_individu_terbaik / 2), 1).*(jumlah_gen - 2)) +1;
point_crossover = point_crossover.*(rand(floor(jumlah_individu_terbaik / 2), 1) < probabilitas_crossover);
for j = 1: length(point_crossover)
	if point_crossover(j)
		chromosome_baru(2*j-1:2*j,:)=[chromosome_baru(2*j-1:2*j,1:point_crossover(j)),...
			flipud(chromosome_baru(2*j-1:2*j,point_crossover(j)+1:jumlah_gen))];
	end
end

% Mutasi
mutasi = find(rand(jumlah_individu_terbaik,jumlah_gen) < probabilitas_mutasi);
chromosome_baru(mutasi)=round(rand(length(mutasi),1));

% Nilai Fitness Baru
xy_baru = zeros(jumlah_individu_terbaik, jumlah_variabel);
for individu = 1:jumlah_variabel
	xy_baru(:,individu) = chromosome_baru(:,1+jumlah_bit_variabel*(individu-1):jumlah_bit_variabel*individu)*[2.^(jumlah_bit_variabel-1:-1:0)]';
	xy_baru(:,individu) = xymin(individu)+(xymax(individu)-xymin(individu))*(xy_baru(:,individu)+1)./(2^jumlah_bit_variabel+1);
end

fx_baru=evalFx(fungsi,xy_baru(:,1),xy_baru(:,2));


% Populasi Baru
if jumlah_individu - jumlah_individu_terbaik
  [ans,Index] = sort(nilai_fitness);
  chromosome = [chromosome(Index(jumlah_individu_terbaik + 1:jumlah_individu), :); chromosome_baru];
  xy = [xy(Index(jumlah_individu_terbaik+1:jumlah_individu),:);xy_baru];
  fx = [fx(Index(jumlah_individu_terbaik+1:jumlah_individu));fx_baru];
else
  chromosome = chromosome_baru;
  x = x_baru;
  fx = fx_baru;
end

% Grafik
mesh(x,y,z)
axis([-3 3 -3 3 0 6])
hold;
contour(x,y,z,20,'k')
hold on;

scatter3(xy(:,1),xy(:,2),fx+4.08,40,'r','filled')
plot(xy(:,1),xy(:,2),'k.','markersize',23)
title(['Generasi # ', num2str(i)]);
xlabel('Parameter "x"');
ylabel('Parameter "y"');
zlabel('Nilai Fitness');
pause(0.2)
 hold;

 best(1+i)=max(fx);
 ave(1+i)=mean(fx);
end

figure('name','Grafik');
plot(0 : jumlah_generasi, best, 0 : jumlah_generasi, ave);
legend('Best', 'Average');
title(['Pc = ', num2str(probabilitas_crossover),', Pm = ', num2str(probabilitas_mutasi)]);
xlabel('Generasi');
ylabel('Nilai Fitness')

function z=evalFx(fx,x,y)
z=eval(fx);