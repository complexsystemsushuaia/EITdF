mdl = mk_thorax_model('male',electrode_positions,[10,10,2],20);

[stim,msel] = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'},1);
opt.imgsz = [32 32];
opt.distr = 3; % non-random, uniform
opt.Nsim = 500; % 500 hundred targets to train on, seems enough
opt.target_size = 0.03; %small targets
opt.target_offset = 0;
opt.noise_figure = .5; % this is key!
opt.square_pixels = 1;
opt.noise_covar = noise_covar;

imdl=mk_GREIT_model(img2, 0.25, [], opt);

imdl.fwd_model.meas_select = msel;
