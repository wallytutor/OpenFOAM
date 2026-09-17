classdef GrafSimulation < handle
    methods (Access = public)
        function [obj] = GrafSimulation(T, P, Y, tout, nsteps, saveas)
            [t_sol, Y_sol] = obj.integrate(T, P, Y, tout, nsteps);

            obj.t_sol = t_sol;
            obj.Y_sol = Y_sol;

            obj.graf_plot(t_sol, 100 * Y_sol, saveas);
        end % function
    end % methods

    methods (Access = private, Static = true)
        function k = arrhenius_rate(T, A, Ea)
            k = A * exp(-Ea / (8.31446261815324 * T));
        end % function arrhenius_rate

        function M = mean_molecular_mass(Y, W)
            M = 1.0 / sum(Y ./ W);
        end % function mean_molecular_mass

        function C = ideal_gas_concentration(rho, Y, W)
            C = rho * Y ./ W;
        end % function ideal_gas_concentration

        function rho = ideal_gas_density(T, P, Y, W)
            M = GrafSimulation.mean_molecular_mass(Y, W);
            rho = (P * M) / (8314.46261815324 * T);
        end % function ideal_gas_density

        function graf_plot(t, Y, saveas)
            % Display solution of kinetics integration in a standard way.
            h = figure();

            subplot(2,2,1);
            plot(t, Y(:,1), "linewidth", 4);
            grid();
            set(gca, 'GridLineStyle', ':');
            xlabel("Time (s)");
            ylabel("Mass percentage");
            l = legend("C_2H_2  ");
            set(l, 'location', 'northeast');

            subplot(2,2,2);
            plot(t, Y(:,[2, 7]), "linewidth", 4);
            grid();
            set(gca, 'GridLineStyle', ':');
            xlabel("Time (s)");
            ylabel("Mass percentage");
            l = legend({"H_2  ", "C_s  "});
            set(l, 'location', 'northwest');

            subplot(2,2,3);
            plot(t, Y(:,[4, 3]), "linewidth", 4);
            grid();
            set(gca, 'GridLineStyle', ':');
            xlabel("Time (s)");
            ylabel("Mass percentage");
            l = legend({"C_2H_4  ", "CH_4 "});
            set(l, 'location', 'northwest');

            subplot(2,2,4);
            plot(t, Y(:,5:6), "linewidth", 4);
            grid();
            set(gca, 'GridLineStyle', ':');
            xlabel("Time (s)");
            ylabel("Mass percentage");
            l = legend({"C_4H_4  ", "C_6H_6  "});
            set(l, 'location', 'northwest');

            print(saveas, '-dpng', '-r300');
        end % function graf_plot
    end % methods

    methods (Access = private, Static = false)
        function [Ydot] = kinetics(self, Y, t, T, P)
            % Compute mass fraction rates of production species.
            % Molar masses C2H2, H2, CH4, C2H4, C4H4, C6H6, Cs, N2.
            W = [26; 2; 16; 28; 52; 78; 12; 28];

            % Matrix of stoichiometric coefficients for the reactions above.
            nu = [
                -1, +1, -1, +1, -2, +2, -1, -1, +0; % 1 C2H2
                -1, +1, -3, +3, +0, +0, +0, +1, +3; % 2 H2
                +0, +0, +2, -2, +0, +0, +0, +0, +0; % 3 CH4
                +1, -1, +0, +0, +0, +0, +0, +0, +0; % 4 C2H4
                +0, +0, +0, +0, +1, -1, -1, +0, +0; % 5 C4H4
                +0, +0, +0, +0, +0, +0, +1, +0, -1; % 6 C6H6
                +0, +0, +0, +0, +0, +0, +0, +2, +6; % 7 Cs
                +0, +0, +0, +0, +0, +0, +0, +0, +0; % 8 N2
            ];

            % Compute gas density.
            rho = 1000.0 * self.ideal_gas_density(T, P, Y, W);

            % Convert mass fractions to concentrations.
            X = self.ideal_gas_concentration(rho, Y, W);

            % Rate constants.
            rate_const = [
                self.arrhenius_rate(T, 4.4e+03, 1.030e+05);
                self.arrhenius_rate(T, 3.8e+07, 2.000e+05);
                self.arrhenius_rate(T, 1.4e+05, 1.500e+05);
                self.arrhenius_rate(T, 8.6e+06, 1.950e+05);
                self.arrhenius_rate(T, 1.2e+05, 1.207e+05);
                self.arrhenius_rate(T, 1.0e+15, 3.352e+05);
                self.arrhenius_rate(T, 1.8e+03, 6.450e+04);
                self.arrhenius_rate(T, 5.5e+06, 1.650e+05);
                self.arrhenius_rate(T, 1.0e+03, 7.500e+04);
            ];

            % Reaction rates in molar units.
            rates = [
                X(1) * X(2)^0.36;
                X(4)^0.50;
                X(1)^0.35 * X(2)^0.22;
                X(3)^0.21;
                X(1)^1.60;
                X(5)^0.75;
                X(1)^1.30 * X(5)^0.60;
                X(1)^1.90 / (1.0 + 18.0 * X(2));
                X(6)^0.75 / (1.0 + 22.0 * X(2));
            ];

            % Rate of molar production of species.
            omega = nu * (rate_const .* rates);

            % Get rate in mass fraction units.
            Wdot = omega .* W / rho;

            % Solving all species at once.
            Ydot = Wdot;
        end % function

        function [ts, Ys] = integrate(self, T, P, Y, tout, nsteps)
            % Integrate kinetics over time and display results.

            % Create anonymous function.
            f =  @(Y, t) self.kinetics(Y, t, T, P);

            ts = linspace(0, tout, nsteps);
            lsode_options('integration method', 'stiff');
            lsode_options('absolute tolerance', 1.0e-15);
            lsode_options('relative tolerance', 1.0e-06);
            lsode_options('minimum step size', 1.0e-12);
            lsode_options('maximum step size', 1.0e-01);
            lsode_options('initial step size', 1.0e-12);
            lsode_options('maximum order', 5);

            [Ys, istate, msg] = lsode(f, Y, ts);
            d = all(abs(sum(Ys') - 1) < 1.0e-06);

            printf("Solution report\n---------------\n")
            printf("- mass conservation: %d\n", d);
            printf("- solver state:      %d\n", istate);
            printf("- status message:    %s\n", msg);
        end
    end % methods

    properties (SetAccess = private, GetAccess = public)
        % Array of time-points corresponding to solution.
        t_sol;

        % Species mass fractions table over time.
        Y_sol;
    end % properties
end % classdef
