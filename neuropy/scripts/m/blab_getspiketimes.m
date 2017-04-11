[animal, series_num, exp_num] = deal('PVCre_0113', 1, 11);

mouse = data.Mice(sprintf('mouse_id="%s"', animal));
series = data.Series(sprintf('series_num=%d', series_num)) & mouse;
exp = data.Experiments(sprintf('exp_num=%d', exp_num)) & series;
units = data.Units & exp;
uids = units.fetchn('unit_id')'; % need to transpose from column vector to row vector for loop

spikes = {};
for uid = uids % 1-based unit IDs
    unit = data.Units(sprintf('unit_id=%d', uid)) & exp;
    ts = unit.unitGetRasters([0], [inf], [0 0]);
    spikes{uid} = ts;
end

save('PVCre_0113_s01_e11_spikes', 'spikes')
