% Optimal Quaternion of Accelerometer-Magnetometer combination (OQAM)
% Proposed by Jin Wu
% e-mail: jin_wu_uestc@hotmail.com


function q = OQAM(Ab, Mb, Mr, w)

    ax = Ab(1);       ay = Ab(2);       az = Ab(3);
    mx = Mb(1);       my = Mb(2);       mz = Mb(3);

    mN = Mr(1);
    mD = Mr(3);
    
    alpha = ax * mx + ay * my + az * mz;
    t2 = sqrt(1 - alpha * alpha);
    t1 = sqrt(1 - 2 * w - 2 * alpha * mD * (w - 1) * w - 2 * t2 * mN * (w - 1) * w + 2 * w^2);

    q = [
           my * (alpha * mN + ax * t1 + mD * t2 * (-1 + w) - alpha * mN * w) + ...
           ay * (-(mx * t1) + mN * (-1 + w) - t2 * w);
         
         - (t1 * t2) + mD * mz * t2 * (1 - w) + alpha * mD * mx * (-1 + w) + alpha * mN * mz * (-1 + w) + ...
           mN * mx * t2 * (-1 + w) - mx * w + ax * (mD - mz * t1 + alpha * w - mD * w) + ...
           az * (mN + mx * t1 - mN * w + t2 * w);

           my * (az * t1 + alpha * mD * (-1 + w) + mN * t2 * (-1 + w) - w) + ...
           ay * (mD - mz * t1 + alpha * w - mD * w);

           alpha * mN * mx * (1 - w) + ax * mN * (-1 + w) + alpha * mD * mz * (-1 + w) + ...
           mD * mx * t2 * (-1 + w) + mN * mz * t2 * (-1 + w) - mz * w - ax * t2 * w + ...
           az * (mD + alpha * w - mD * w);
    ];

    q = q ./ norm(q);
end