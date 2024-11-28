%% For loop demonstration

demo_elements = ["C","H","O","N","Se"];
demo_mass = [12,1,16,14,80];

% Compound CH4
% Manual method
demo_composition = [1,4,0,0,1];

demo_mass_total = demo_mass(1)*demo_composition(1)+demo_mass(2)*demo_composition(2)+demo_mass(3)*demo_composition(3)+demo_mass(4)*demo_composition(4); %(12*1)+(1*4)+(0*16)+(0*14)=16


% Sequential method
demo_mass_total2 = 0;
demo_mass_total2 = demo_mass_total2+demo_mass(1)*demo_composition(1); %  0 + (12*1) = 12
demo_mass_total2 = demo_mass_total2+demo_mass(2)*demo_composition(2); % 12 + (1*4) = 16
demo_mass_total2 = demo_mass_total2+demo_mass(3)*demo_composition(3); % 16 + (0*16) = 16
demo_mass_total2 = demo_mass_total2+demo_mass(4)*demo_composition(4); % 16 + (0*14) = 16

% Programming method
demo_mass_total3 = 0;
for i=1:size(demo_composition,2)
    demo_mass_total3 = demo_mass_total3+demo_mass(i)*demo_composition(i); %    
end
