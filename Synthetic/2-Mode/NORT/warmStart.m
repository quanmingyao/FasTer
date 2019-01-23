function [ R ] = warmStart(U1, U0, maxR)

R = [U1, U0];
[R, ~] = qr(R, 0);
R = R(:,1:min(size(R,2), maxR));

end