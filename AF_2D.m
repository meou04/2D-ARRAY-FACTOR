classdef AF_2D < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        GridLayout                    matlab.ui.container.GridLayout
        LeftPanel                     matlab.ui.container.Panel
        FrequencyMHzEditField_2Label  matlab.ui.control.Label
        FrequencyMHzEditField         matlab.ui.control.NumericEditField
        AntennaelementDropDownLabel   matlab.ui.control.Label
        AntennaelementDropDown        matlab.ui.control.DropDown
        ArrayConfigurationLabel       matlab.ui.control.Label
        XLabel                        matlab.ui.control.Label
        AxisLabel                     matlab.ui.control.Label
        YLabel                        matlab.ui.control.Label
        ElementsLabel                 matlab.ui.control.Label
        dLabel                        matlab.ui.control.Label
        PhaseshiftLabel               matlab.ui.control.Label
        BeamsteeringLabel             matlab.ui.control.Label
        FeedingcoefficientsDropDownLabel  matlab.ui.control.Label
        FeedingcoefficientsDropDown   matlab.ui.control.DropDown
        PlotButton                    matlab.ui.control.Button
        DirectivitydBiEditFieldLabel  matlab.ui.control.Label
        DirectivitydBiEditField       matlab.ui.control.NumericEditField
        PhasetypeDropDownLabel        matlab.ui.control.Label
        PhasetypeDropDown             matlab.ui.control.DropDown
        Label                         matlab.ui.control.Label
        EditField                     matlab.ui.control.NumericEditField
        EditField_2Label              matlab.ui.control.Label
        EditField_2                   matlab.ui.control.NumericEditField
        Elements_X                    matlab.ui.control.NumericEditField
        Elements_Y                    matlab.ui.control.NumericEditField
        DLambda_X                     matlab.ui.control.NumericEditField
        DLambda_Y                     matlab.ui.control.NumericEditField
        Phasestep_X                   matlab.ui.control.NumericEditField
        Phasestep_Y                   matlab.ui.control.NumericEditField
        RightPanel                    matlab.ui.container.Panel
        UIAxes                        matlab.ui.control.UIAxes
        DPlannarArraySimulatorLabel   matlab.ui.control.Label
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
    end

    
    methods (Access = private)
        
        function updateAntennaPhase(app)
            phase_type = app.PhasetypeDropDown.Value;
            frequency = app.FrequencyMHzEditField.Value;
            c = 3e8;
            wavelength = c / (frequency * 1e6);
            wavelength_mm = wavelength * 1000;
            k = 2*pi / wavelength_mm;
            dx = app.DLambda_X.Value*wavelength_mm;
            dy = app.DLambda_Y.Value*wavelength_mm;
            
            switch phase_type
                case 'Manual'
                    DP_x = deg2rad(app.Phasestep_X.Value);
                    DP_y = deg2rad(app.Phasestep_Y.Value);
                    u0 = -DP_x / (k * dx);
                    v0 = -DP_y / (k * dy);
                    
                    sin_theta_sq = u0^2 + v0^2;
                    theta_calc = rad2deg(asin(sqrt(sin_theta_sq)));
                    phi_calc = rad2deg(atan2(v0, u0));
                    
                    app.EditField.Value = theta_calc;
                    app.EditField_2.Value = phi_calc;
                    
                case 'Steering'
                    theta_rad = deg2rad(app.EditField.Value);
                    phi_rad   = deg2rad(app.EditField_2.Value);
                    u0 = sin(theta_rad) * cos(phi_rad);
                    v0 = sin(theta_rad) * sin(phi_rad);
                    
                    DP_x_disp = rad2deg(-k * dx * u0);
                    DP_y_disp = rad2deg(-k * dy * v0);
                    
                    app.Phasestep_X.Value = DP_x_disp;
                    app.Phasestep_Y.Value = DP_y_disp;
            end
        end
        
        
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.PhasetypeDropDownValueChanged([]);
        end

        % Button pushed function: PlotButton
        function PlotButtonPushed(app, event)
            f_MHz = app.FrequencyMHzEditField.Value;
            lambda_mm = 3e8 / (f_MHz * 1e6) * 1000; 
            k = 2*pi / lambda_mm;
            Nx = app.Elements_X.Value;
            Ny = app.Elements_Y.Value;
            dx = app.DLambda_X.Value*lambda_mm;
            dy = app.DLambda_Y.Value*lambda_mm;
            scan_theta = app.EditField.Value;
            scan_phi = app.EditField_2.Value;
            theta_rad = deg2rad(scan_theta);
            phi_rad   = deg2rad(scan_phi);
            u0 = sin(theta_rad) * cos(phi_rad);
            v0 = sin(theta_rad) * sin(phi_rad);
            
            feeding_type = app.FeedingcoefficientsDropDown.Value;
            switch feeding_type
                case 'Uniform'
                    wx = ones(Nx, 1);
                    wy = ones(Ny, 1);
                case 'Binomial'
                    if Nx > 1, wx = diag(rot90(pascal(Nx)))'; else, wx = 1; end
                    if Ny > 1, wy = diag(rot90(pascal(Ny)))'; else, wy = 1; end
                case 'Triangular'
                    wx = triang(Nx);
                    wy = triang(Ny);
            end
            
            theta_plot = linspace(0, 90, 91);
            phi_plot   = linspace(0, 360, 361);
            [THETA, PHI] = meshgrid(deg2rad(theta_plot), deg2rad(phi_plot));
            U = sin(THETA) .* cos(PHI);
            V = sin(THETA) .* sin(PHI);
            
            antenna_element_type = app.AntennaelementDropDown.Value;
            switch antenna_element_type
                case 'Isotropic'
                    element_pattern = ones(size(THETA));
                case 'lambda/2_dipole'
                    element_pattern = abs((cos(pi/2 * sin(THETA)) + eps) ./ (cos(THETA) + eps));
                    element_pattern(isnan(element_pattern)) = 0;
            end
                   
            AF_x = zeros(size(THETA));
            for m = 0:Nx-1
                phase_x_steer = -k * m * dx * u0;
                AF_x = AF_x + wx(m+1) * exp(1j * (k * m * dx * U + phase_x_steer));
            end
            AF_y = zeros(size(THETA));
            for n = 0:Ny-1
                phase_y_steer = -k * n * dy * v0;
                AF_y = AF_y + wy(n+1) * exp(1j * (k * n * dy * V + phase_y_steer));
            end
            AF = AF_x .* AF_y;
            AF_mag = abs(AF);
            AF_norm = AF_mag / max(AF_mag(:));
            
            Total_Pattern = AF_norm .* element_pattern;
            Total_Pattern_norm = Total_Pattern / max(Total_Pattern(:));
            
            Total_Pattern_dB = 20*log10(Total_Pattern_norm + 1e-6);
            Total_Pattern_dB(Total_Pattern_dB < -40) = -40;
            
            R = Total_Pattern_norm; 
            X = R .* sin(THETA) .* cos(PHI);
            Y = R .* sin(THETA) .* sin(PHI);
            Z = R .* cos(THETA);
            
            surf(app.UIAxes,X, Y, Z, Total_Pattern_dB);
            shading(app.UIAxes,'interp');
            colormap(app.UIAxes, 'jet');
            cb = colorbar(app.UIAxes);
            caxis(app.UIAxes,[-40 0]);
            xlabel(app.UIAxes,'X'); ylabel(app.UIAxes,'Y'); zlabel(app.UIAxes,'Z');
            axis(app.UIAxes, 'equal');
            grid(app.UIAxes, 'on');
            view(app.UIAxes, 45, 30);
            
            U = Total_Pattern_norm.^2;
            U_max = max(U(:));
            d_theta = deg2rad(theta_plot(2) - theta_plot(1)); % Bÿÿc nhÿy theta (rad)
            d_phi   = deg2rad(phi_plot(2) - phi_plot(1));     % Bÿÿc nhÿy phi (rad)
            P_rad = sum(sum(U .* sin(THETA))) * d_theta * d_phi;
            Directivity_linear = 4 * pi * U_max / P_rad;
            Directivity_dB = 10 * log10(Directivity_linear); 
            app.DirectivitydBiEditField.Value = Directivity_dB;
        end

        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.UIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 2x1 grid
                app.GridLayout.RowHeight = {741, 741};
                app.GridLayout.ColumnWidth = {'1x'};
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 1;
            else
                % Change to a 1x2 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {343, '1x'};
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 2;
            end
        end

        % Value changed function: PhasetypeDropDown
        function PhasetypeDropDownValueChanged(app, event)
            phase_type = app.PhasetypeDropDown.Value;
            switch phase_type
                case 'Manual'
                    app.Phasestep_X.Enable = 'on';
                    app.Phasestep_Y.Enable = 'on';
                    app.EditField.Enable = 'off';
                    app.EditField_2.Enable = 'off';
                case 'Steering'
                    app.Phasestep_X.Enable = 'off';
                    app.Phasestep_Y.Enable = 'off';
                    app.EditField.Enable = 'on';
                    app.EditField_2.Enable = 'on';
            end
            updateAntennaPhase(app);
        end

        % Value changed function: Phasestep_X
        function Phasestep_XValueChanged(app, event)
            updateAntennaPhase(app);
        end

        % Value changed function: Phasestep_Y
        function Phasestep_YValueChanged(app, event)
            updateAntennaPhase(app);
        end

        % Value changed function: EditField
        function EditFieldValueChanged(app, event)
            updateAntennaPhase(app);
        end

        % Value changed function: EditField_2
        function EditField_2ValueChanged(app, event)
            updateAntennaPhase(app);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Position = [100 100 1302 741];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {343, '1x'};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            % Create FrequencyMHzEditField_2Label
            app.FrequencyMHzEditField_2Label = uilabel(app.LeftPanel);
            app.FrequencyMHzEditField_2Label.HorizontalAlignment = 'right';
            app.FrequencyMHzEditField_2Label.FontSize = 14;
            app.FrequencyMHzEditField_2Label.Position = [8 350 114 22];
            app.FrequencyMHzEditField_2Label.Text = 'Frequency(MHz):';

            % Create FrequencyMHzEditField
            app.FrequencyMHzEditField = uieditfield(app.LeftPanel, 'numeric');
            app.FrequencyMHzEditField.FontSize = 14;
            app.FrequencyMHzEditField.Position = [207 350 128 22];
            app.FrequencyMHzEditField.Value = 1000;

            % Create AntennaelementDropDownLabel
            app.AntennaelementDropDownLabel = uilabel(app.LeftPanel);
            app.AntennaelementDropDownLabel.HorizontalAlignment = 'right';
            app.AntennaelementDropDownLabel.FontSize = 14;
            app.AntennaelementDropDownLabel.Position = [6 663 115 22];
            app.AntennaelementDropDownLabel.Text = 'Antenna element:';

            % Create AntennaelementDropDown
            app.AntennaelementDropDown = uidropdown(app.LeftPanel);
            app.AntennaelementDropDown.Items = {'Isotropic', 'lambda/2_dipole'};
            app.AntennaelementDropDown.FontSize = 14;
            app.AntennaelementDropDown.Position = [164 663 169 22];
            app.AntennaelementDropDown.Value = 'lambda/2_dipole';

            % Create ArrayConfigurationLabel
            app.ArrayConfigurationLabel = uilabel(app.LeftPanel);
            app.ArrayConfigurationLabel.FontSize = 14;
            app.ArrayConfigurationLabel.Position = [11 632 126 22];
            app.ArrayConfigurationLabel.Text = 'Array Configuration';

            % Create XLabel
            app.XLabel = uilabel(app.LeftPanel);
            app.XLabel.FontSize = 14;
            app.XLabel.Position = [38 568 13 22];
            app.XLabel.Text = 'X';

            % Create AxisLabel
            app.AxisLabel = uilabel(app.LeftPanel);
            app.AxisLabel.FontSize = 14;
            app.AxisLabel.Position = [28 601 32 22];
            app.AxisLabel.Text = 'Axis';

            % Create YLabel
            app.YLabel = uilabel(app.LeftPanel);
            app.YLabel.FontSize = 14;
            app.YLabel.Position = [37 536 14 22];
            app.YLabel.Text = 'Y';

            % Create ElementsLabel
            app.ElementsLabel = uilabel(app.LeftPanel);
            app.ElementsLabel.FontSize = 14;
            app.ElementsLabel.Position = [85 601 64 22];
            app.ElementsLabel.Text = 'Elements';

            % Create dLabel
            app.dLabel = uilabel(app.LeftPanel);
            app.dLabel.FontSize = 14;
            app.dLabel.Position = [177 601 25 22];
            app.dLabel.Text = 'd/ÿ';

            % Create PhaseshiftLabel
            app.PhaseshiftLabel = uilabel(app.LeftPanel);
            app.PhaseshiftLabel.HorizontalAlignment = 'right';
            app.PhaseshiftLabel.FontSize = 14;
            app.PhaseshiftLabel.Position = [239 601 97 22];
            app.PhaseshiftLabel.Text = 'Phase shift (º):';

            % Create BeamsteeringLabel
            app.BeamsteeringLabel = uilabel(app.LeftPanel);
            app.BeamsteeringLabel.FontSize = 14;
            app.BeamsteeringLabel.Position = [13 454 100 22];
            app.BeamsteeringLabel.Text = 'Beam steering:';

            % Create FeedingcoefficientsDropDownLabel
            app.FeedingcoefficientsDropDownLabel = uilabel(app.LeftPanel);
            app.FeedingcoefficientsDropDownLabel.VerticalAlignment = 'top';
            app.FeedingcoefficientsDropDownLabel.FontSize = 14;
            app.FeedingcoefficientsDropDownLabel.Position = [9 385 141 22];
            app.FeedingcoefficientsDropDownLabel.Text = '  Feeding coefficients:';

            % Create FeedingcoefficientsDropDown
            app.FeedingcoefficientsDropDown = uidropdown(app.LeftPanel);
            app.FeedingcoefficientsDropDown.Items = {'Uniform', 'Binomial', 'Triangular'};
            app.FeedingcoefficientsDropDown.FontSize = 14;
            app.FeedingcoefficientsDropDown.Position = [167 385 169 22];
            app.FeedingcoefficientsDropDown.Value = 'Uniform';

            % Create PlotButton
            app.PlotButton = uibutton(app.LeftPanel, 'push');
            app.PlotButton.ButtonPushedFcn = createCallbackFcn(app, @PlotButtonPushed, true);
            app.PlotButton.FontSize = 15;
            app.PlotButton.Position = [38 202 262 48];
            app.PlotButton.Text = 'Plot';

            % Create DirectivitydBiEditFieldLabel
            app.DirectivitydBiEditFieldLabel = uilabel(app.LeftPanel);
            app.DirectivitydBiEditFieldLabel.HorizontalAlignment = 'right';
            app.DirectivitydBiEditFieldLabel.FontSize = 14;
            app.DirectivitydBiEditFieldLabel.Position = [9 150 104 22];
            app.DirectivitydBiEditFieldLabel.Text = 'Directivity (dBi):';

            % Create DirectivitydBiEditField
            app.DirectivitydBiEditField = uieditfield(app.LeftPanel, 'numeric');
            app.DirectivitydBiEditField.FontSize = 14;
            app.DirectivitydBiEditField.Enable = 'off';
            app.DirectivitydBiEditField.Position = [208 150 128 22];

            % Create PhasetypeDropDownLabel
            app.PhasetypeDropDownLabel = uilabel(app.LeftPanel);
            app.PhasetypeDropDownLabel.HorizontalAlignment = 'right';
            app.PhasetypeDropDownLabel.FontSize = 14;
            app.PhasetypeDropDownLabel.Position = [6 497 80 22];
            app.PhasetypeDropDownLabel.Text = ' Phase type:';

            % Create PhasetypeDropDown
            app.PhasetypeDropDown = uidropdown(app.LeftPanel);
            app.PhasetypeDropDown.Items = {'Manual', 'Steering'};
            app.PhasetypeDropDown.ValueChangedFcn = createCallbackFcn(app, @PhasetypeDropDownValueChanged, true);
            app.PhasetypeDropDown.FontSize = 14;
            app.PhasetypeDropDown.Position = [164 497 169 22];
            app.PhasetypeDropDown.Value = 'Steering';

            % Create Label
            app.Label = uilabel(app.LeftPanel);
            app.Label.HorizontalAlignment = 'right';
            app.Label.FontSize = 14;
            app.Label.Position = [136 454 35 22];
            app.Label.Text = ' ÿ (º):';

            % Create EditField
            app.EditField = uieditfield(app.LeftPanel, 'numeric');
            app.EditField.ValueChangedFcn = createCallbackFcn(app, @EditFieldValueChanged, true);
            app.EditField.Position = [197 454 100 22];

            % Create EditField_2Label
            app.EditField_2Label = uilabel(app.LeftPanel);
            app.EditField_2Label.HorizontalAlignment = 'right';
            app.EditField_2Label.FontSize = 14;
            app.EditField_2Label.Position = [135 422 36 22];
            app.EditField_2Label.Text = {'ÿ (º):'; ''};

            % Create EditField_2
            app.EditField_2 = uieditfield(app.LeftPanel, 'numeric');
            app.EditField_2.ValueChangedFcn = createCallbackFcn(app, @EditField_2ValueChanged, true);
            app.EditField_2.Position = [197 422 100 22];

            % Create Elements_X
            app.Elements_X = uieditfield(app.LeftPanel, 'numeric');
            app.Elements_X.Position = [85 568 51 22];
            app.Elements_X.Value = 5;

            % Create Elements_Y
            app.Elements_Y = uieditfield(app.LeftPanel, 'numeric');
            app.Elements_Y.Position = [85 536 51 22];
            app.Elements_Y.Value = 5;

            % Create DLambda_X
            app.DLambda_X = uieditfield(app.LeftPanel, 'numeric');
            app.DLambda_X.Position = [164 568 51 22];
            app.DLambda_X.Value = 0.5;

            % Create DLambda_Y
            app.DLambda_Y = uieditfield(app.LeftPanel, 'numeric');
            app.DLambda_Y.Position = [164 536 51 22];
            app.DLambda_Y.Value = 0.5;

            % Create Phasestep_X
            app.Phasestep_X = uieditfield(app.LeftPanel, 'numeric');
            app.Phasestep_X.ValueChangedFcn = createCallbackFcn(app, @Phasestep_XValueChanged, true);
            app.Phasestep_X.Position = [247 568 51 22];

            % Create Phasestep_Y
            app.Phasestep_Y = uieditfield(app.LeftPanel, 'numeric');
            app.Phasestep_Y.ValueChangedFcn = createCallbackFcn(app, @Phasestep_YValueChanged, true);
            app.Phasestep_Y.Position = [247 536 51 22];

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;
            app.RightPanel.FontSize = 14;

            % Create UIAxes
            app.UIAxes = uiaxes(app.RightPanel);
            title(app.UIAxes, '2D Radiation Pattern')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            app.UIAxes.FontSize = 14;
            app.UIAxes.Position = [15 55 934 590];

            % Create DPlannarArraySimulatorLabel
            app.DPlannarArraySimulatorLabel = uilabel(app.RightPanel);
            app.DPlannarArraySimulatorLabel.HorizontalAlignment = 'center';
            app.DPlannarArraySimulatorLabel.FontName = 'Tahoma';
            app.DPlannarArraySimulatorLabel.FontSize = 25;
            app.DPlannarArraySimulatorLabel.FontWeight = 'bold';
            app.DPlannarArraySimulatorLabel.Position = [39 699 355 33];
            app.DPlannarArraySimulatorLabel.Text = {'2D Plannar Array Simulator'; ''};

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = AF_2D

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end