model_name = "male_16el_adjacent_"

task = 1; %0: simmulation. 1:mining.

noise_covar   = 1

normalize_flag  = 0; % 0: No normalizado, 1: normalizado

%
% ELECTRODOS, ESTIMULACIÓN Y MEDICIONES
%
nelec  = 16; % Número de electrodos
nrings = 1;  % Número de anillos de electrodos
electrode_shape = [10,10,2]; % Forma del electrodo (lado x lado x espesor en mm)
electrode_positions = [ ...  % Coordenadas cartesianas, x,y,z, del modelo de torax "male"
   -34.4564,   153.6012,   210.0000; ...
  -103.7541,   150.4564,   210.0000; ...
  -171.2636,   137.0363,   210.0000; ...
  -188.9420,    78.2578,   210.0000; ...
  -176.3317,    10.6117,   210.0000; ...
  -161.7638,   -55.8329,   210.0000; ...
  -112.8155,  -103.2164,   210.0000; ...
   -44.4994,  -108.8109,   210.0000; ...
    24.0444,  -106.0791,   210.0000; ...
    92.7143,  -109.3225,   210.0000; ...
   148.8393,   -70.6724,   210.0000; ...
   169.0922,    -5.0908,   210.0000; ...
   167.8353,    64.0442,   210.0000; ...
   168.4765,   131.5265,   210.0000; ...
   104.1889,   150.9233,   210.0000; ...
    34.8602,   153.2953,   210.0000];
amp_level       = 1; %  Corriente de estimulación en miliAmperes
inj  = 'ad';  % Adjacent injection, equivalent to [0,1]
meas = 'ad';  % Adjacent measurement, equivalent to [0,1]
meas_current    = 1; % 1: meas_current. Else: no_meas_current
if (meas_current == 1)
  meas_flag = 'meas_current'
else
  meas_flag = 'no_meas_current'
end

%
% MALLADO
%
max_elem_size   = 6; % Tamaño máximo del elemento del mallado


%
% CONDUCTIVIDADES
%
bkground_cond = 0.33; % En S/m, sacados del paper Consensus
lung_cond_min = 0.042; % En S/m, sacados del paper Consensus
lung_cond_max = 0.11; % En S/m, sacados del paper Consensus
hearth_cond = 0.48; % En S/m, sacados del paper Consensus


%
% PARAMETROS SIMULACIÓN
%
lung_freq = 0.6; % En Hz
sampling_rate = 12; % En Hz
Duration = 4; % En segundos.
Phase = 1.63 * pi; % Fase de la sinusoide (la obtuve aleatoria, pero la dejo fija para repetibilidad)


%
% MODELO FORWARD
%


%
% MODELO INVERSO
%

% Lista de paths modificables por el usuario:
datadir = '/home/daniel/ORDENAR/UNTDF/SistemasComplejos/Bioimpedancia/GREITdF/sandbox';
eidors  = '/home/daniel/ORDENAR/UNTDF/SistemasComplejos/Bioimpedancia/GREITdF/eidors_mining_version/eidors/'
sandbox = datadir;

%% Paso 1: Inicializar

cd(eidors);
run startup.m;
cd(sandbox);

nelec=size(electrode_positions,1); % Numero de electrodos

%% Paso 2: Hacer el mallado y generar el patron de medición
% Obs: Si se guarda el resultado, se puede saltear este paso

mdl = mk_thorax_model('male',electrode_positions,electrode_shape, max_elem_size); 

mdl.stimulation =  mk_stim_patterns(nelec,1,[0,1],[0,1],{'meas_current'}, amp_level);
mdl = mdl_normalize(mdl, normalize_flag); % Normalizar el modelo o no


% Setting contact impedances
for i = 1:16
   mdl.electrode(i).z_contact = i*1000000000;
end    
mdl.stimulation =  mk_stim_patterns(nelec,1,[0,1],[0,1],{meas_flag}, amp_level);
mdl = mdl_normalize(mdl, normalize_flag); % Normalizar el modelo o no

%% Paso 2.5: generar unas gráficas para ver el modelo.
% Generada con conductividades promedio

% load modelo_06.mat mdl % Opcional: cargar el modelo si fue guardado

clear opt % Borrar el vector de opciones:
opt.bkgnd_val = bkground_cond;
opt.lung_val = (lung_cond_min + lung_cond_max) / 2;
opt.heart_val = hearth_cond;
% Correr un toque a la izquierda los órganos
opt.left_lung_center = [71.648 15.0871 88.062]; % Default: [61.648 5.0871 88.062]
opt.right_lung_center = [-71.47 15.0871 88.062]; % Default: [-81.47 5.0871 88.062]
opt.heart_center = [20.96 -34.607 187.45]; % Default: [10.96 -44.607 187.45]

img = mk_lung_image(mdl,opt); % Generar modelo con conductividades

% Modelo, pulmón y corazón
figure; show_fem(img); % imagen generada por EIDORS

save pulmones.mat img

% Corte transversal a la altura de los electrodos
figure; vis3d_patchFlat(img.fwd_model.nodes,img.fwd_model.elems(:,1:4),img.elem_data,'p(:,3)<210');
colormap('jet'); 
colorbar; material dull; 
%lighting phong; 
camlight; % Para que la imagen quede mas linda

% Pulmones
figure; vis3d(img.fwd_model.nodes,img.fwd_model.elems(img.elem_data<=lung_cond_max,1:4),img.fwd_model.nodes(:,1)); 
colormap('jet'); 
colorbar; 
material dull; 
%lighting phong; 
camlight; % Para que la imagen quede mas linda

disp("**************************************");
disp("**     GENERANDO SINUSOIDE          **");
disp("**************************************");

%% Paso 3: generar la sinusoide:
t=0:1/sampling_rate:Duration; % time samples in seconds
lung_cond_dynamic = (lung_cond_max+lung_cond_min)/2+...
    (lung_cond_max-lung_cond_min)/2*sin(2*pi*lung_freq*t+Phase); % Sinusoide de conductividad
% figure; plot(t,lung_cond_dynamic); % Chequeo

disp("**************************************");
disp("**  GENERANDO MUESTRAS TEMPORALES   **");
disp("**************************************");

%% Paso 4: generar las muestras temporales
if (meas_current == 1)
  size_frame = nelec*nelec;
else
  size_frame = (nelec-3)*nelec;
end
output_signals = zeros (size_frame,length(t)); % Inicializo la variable de medidas
for rr=1:length(t)
    tic; % Empezar cronometro
    opt.lung_val = lung_cond_dynamic(rr); % modificar cond del pulmón.
    img = mk_lung_image(mdl,opt); % Generar las conductividades
    vh = fwd_solve(img); % Resolver el PD
    output_signals(:,rr) = vh.meas; % Guardar las señales simuladas
    t_elapsed = toc; % Finalizar cronometro
    display(['Iteración ' num2str(rr) ' de ' num2str(length(t)) '. Tiempo de cómputo = ' num2str(t_elapsed)]);
end

% Graficar señales
%figure; plot(t,output_signals);

disp("**************************************");
disp("**         AGREGANDO RUIDO          **");
disp("**************************************");

#FRAMES = (add_noise(500,output_signals)).meas;
FRAMES = (add_noise(5000,output_signals)).meas;

% Guardar senales generadas, modelo y parámetros (opcional)
% save simulated_signals_v1.mat output_signals mdl amp_level electrode_positions opt Phase sampling_rate t lung_freq normalize_flag electrode_shape max_elem_size *cond*

disp("**************************************");
disp("**            GUARDANDO             **");
disp("**************************************");

save -6 simuresp_ztest.mat FRAMES
%save simuresp_ztest_rw.mat output_signals

disp("**************************************");
disp("**            LEEESTOOO             **");
disp("**************************************");

exit;

disp("**************************************");
disp("**  MISTERY!!!!!!!!!!!!!!!!!!!!!!   **");
disp("**************************************");
