function GA_1

disp ('Algoritma Genetika : Fungsi dengan 1 variable')
disp (' ')
disp ('Masalah : Mencari nilai maksimum dari fungsi (15*x - x*x) dimana parameter "x" berada pada rang 0 dan 15')

% Inisialisasi Populasi
jumlah_individu = 6;
jumlah_gen = 4;
xmin = 0;
xmax = 15;
jumlah_generasi = 20;
probabilitas_crossover = 0.9;
probabilitas_mutasi = 0.001;


% pupulasi awal
chromosome = round(rand(jumlah_individu, jumlah_gen));

% menghitung nilai x (konversi biner ke desimal)
x = chromosome * (2.^(jumlah_gen-1:-1:0))';

fungsi = '15*x -x.^2';

% Menghitung nilai fx
fx = evalFx(fungsi,x);
best(1) = max(fx);
ave(1) = mean(fx);

% Menampilkan grafik
figure('name', 'Lokasi chromosome')
fplot(fungsi, [xmin,xmax])
hold;
plot(x, fx, 'r. ', 'markersize', 15)
legend(['fungsi : ', fungsi], 'Populasi awal');
title('Hit any key to continue');
xlabel('Parameter "x"');
ylabel('Nilai fitness');
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
[selind,~] = find(chance < cumulative_fitness & chance >= [zeros(1, jumlah_individu_terbaik); cumulative_fitness(1: jumlah_individu-1, :)]);
chromosome_baru = chromosome(selind, :);

% Crossover
point_crossover = round(rand(floor(jumlah_individu_terbaik / 2), 1).*(jumlah_gen - 2)) +1;
point_crossover = point_crossover.*(rand(floor(jumlah_individu_terbaik / 2), 1) < probabilitas_crossover);
for j=1:length(point_crossover)
    if point_crossover(j)   
	chromosome_baru(2*j-1:2*j,:)=[chromosome_baru(2*j-1:2*j,1:point_crossover(j)),...
		flipud(chromosome_baru(2*j-1:2*j,point_crossover(j)+1:jumlah_gen))];
    end
end

% Mutasi
mutasi = find(rand(jumlah_individu_terbaik, jumlah_gen) < probabilitas_mutasi);
chromosome_baru(mutasi) = round(rand(length(mutasi),1));

% Nilai Fitness Baru
x_baru = chromosome_baru * (2.^(jumlah_gen-1: -1: 0))';
x_baru = xmin + (x_baru +1)*(xmax-xmin)/(2^jumlah_gen-1);
fx_baru = evalFx(fungsi,x_baru);

% Populasi Baru
if jumlah_individu - jumlah_individu_terbaik
   [ans,Index]=sort(nilai_fitness);
   chromosome = [chromosome(Index(jumlah_individu_terbaik + 1: jumlah_individu),:); chromosome_baru];
   x = [x(Index(jumlah_individu_terbaik+1: jumlah_individu)); x_baru];
   fx = [fx(Index(jumlah_individu_terbaik+1: jumlah_individu)); fx_baru];
else
   chromosome = chromosome_baru;
   x = x_baru;
   fx = fx_baru;
end

% Menampilkan Grafik
fplot(fungsi, [xmin,xmax])
hold;
plot(x, fx, 'r.', 'markersize', 15)
legend(['Fungsi: ', fungsi], 'Populasi sekarang');
title(['Genarasi # ',num2str(i)]);
xlabel('Parameter "x:');
ylabel('Nilai fitness');
pause(0.2)
hold;

best(1+i) = max(fx);
ave(1+i) = mean(fx);
end

figure('name', 'Grafik');
plot(0 : jumlah_generasi, best, 0 : jumlah_generasi, ave);
legend('Best', 'Avarege');
title(['Pc = ', num2str(probabilitas_crossover),', Pm =',num2str(probabilitas_mutasi)]);
xlabel('Generasi');
ylabel('Nilai Fitness');

function y=evalFx(fx,x)
y=eval(fx);

end
end