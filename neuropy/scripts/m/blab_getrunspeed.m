bm = BallMovement('/Volumes/tudata/ephys/Ntsr1-Cre_0174/174-2-005.nev')

roll = 100 * bm.roll(:)';
pitch = 100 * bm.pitch(:)';
yaw = 100 * bm.yaw(:)';

tspeed = bm.timestamps;
speed = sqrt(roll .^ 2 + pitch .^ 2 + yaw .^ 2);
save('Ntsr1-Cre_0174_s02_e05_runspeed', 'tspeed', 'speed');
