function save_to_csv()
    % rooftop.mat 파일을 불러옴
    clear;
    
    data = load('../data/rooftop.mat');
    
    % Rover 객체의 데이터를 CSV로 저장
    save_rover_data_to_csv(data.Station, '../data/data_station/');
    save_rover_data_to_csv(data.Rover, '../data/data_rover/');
end

function save_rover_data_to_csv(rover, output_dir)
    % output_dir 경로가 없으면 생성
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    % 필요한 데이터 추출 및 CSV 파일로 저장
    save_pr1_to_csv(rover, fullfile(output_dir, 'pr.csv'));
    save_dop1_to_csv(rover, fullfile(output_dir, 'dop.csv'));
    save_carrier1_to_csv(rover, fullfile(output_dir), 'carrier.csv');
    save_sv_pos_to_csv(rover, fullfile(output_dir, 'sv_pos.csv'));
    save_sv_vel_to_csv(rover, fullfile(output_dir, 'sv_vel.csv'));
    save_gps_data_to_csv(rover, fullfile(output_dir, 'time.csv'));
end

function save_carrier1_to_csv(rover, output_dir, filename)
    % carrier 데이터 추출
    carrier = rover.ph1;

    % 테이블로 변환
    T = table(carrier);

    % 파일 전체 경로 생성
    output_filename = fullfile(output_dir, filename);

    % 테이블을 CSV 파일로 저장
    writetable(T, output_filename);
end

function save_pr1_to_csv(rover, output_filename)
    % pr1 데이터 추출
    pr1 = rover.pr1;

    % 테이블로 변환
    T = table(pr1);

    % 테이블을 CSV 파일로 저장
    writetable(T, output_filename);
end

function save_dop1_to_csv(rover, output_filename)
    % dop1 데이터 추출
    dop1 = rover.dop1;

    % 테이블로 변환
    T = table(dop1);

    % 테이블을 CSV 파일로 저장
    writetable(T, output_filename);
end

function save_sv_pos_to_csv(rover, output_filename)
    % SVpos_x, SVpos_y, SVpos_z 데이터 추출
    SVpos_x = rover.SVpos_x;
    SVpos_y = rover.SVpos_y;
    SVpos_z = rover.SVpos_z;

    % 데이터 크기 확인
    [num_samples, num_satellites] = size(SVpos_x);
    
    % 데이터 저장을 위한 테이블 생성
    var_names = cell(1, num_satellites * 3);
    data = zeros(num_samples, num_satellites * 3);
    
    for i = 1:num_satellites
        % 변수 이름 설정
        var_names{(i-1)*3 + 1} = sprintf('SVpos_x_%d', i);
        var_names{(i-1)*3 + 2} = sprintf('SVpos_y_%d', i);
        var_names{(i-1)*3 + 3} = sprintf('SVpos_z_%d', i);
        
        % 데이터 설정
        data(:, (i-1)*3 + 1) = SVpos_x(:, i);
        data(:, (i-1)*3 + 2) = SVpos_y(:, i);
        data(:, (i-1)*3 + 3) = SVpos_z(:, i);
    end
    
    % 테이블로 변환
    T = array2table(data, 'VariableNames', var_names);

    % 테이블을 CSV 파일로 저장
    writetable(T, output_filename);
end

function save_sv_vel_to_csv(rover, output_filename)
    % SVvel_x, SVvel_y, SVvel_z 데이터 추출
    SVvel_x = rover.SVvel_x;
    SVvel_y = rover.SVvel_y;
    SVvel_z = rover.SVvel_z;

    % 데이터 크기 확인
    [num_samples, num_satellites] = size(SVvel_x);
    
    % 데이터 저장을 위한 테이블 생성
    var_names = cell(1, num_satellites * 3);
    data = zeros(num_samples, num_satellites * 3);
    
    for i = 1:num_satellites
        % 변수 이름 설정
        var_names{(i-1)*3 + 1} = sprintf('SVvel_x_%d', i);
        var_names{(i-1)*3 + 2} = sprintf('SVvel_y_%d', i);
        var_names{(i-1)*3 + 3} = sprintf('SVvel_z_%d', i);
        
        % 데이터 설정
        data(:, (i-1)*3 + 1) = SVvel_x(:, i);
        data(:, (i-1)*3 + 2) = SVvel_y(:, i);
        data(:, (i-1)*3 + 3) = SVvel_z(:, i);
    end
    
    % 테이블로 변환
    T = array2table(data, 'VariableNames', var_names);

    % 테이블을 CSV 파일로 저장
    writetable(T, output_filename);
end

function save_gps_data_to_csv(rover, output_filename)
    % time_GPS와 week 데이터 추출
    gps_week = rover.week;
    gps_time = rover.time_GPS;

    % 테이블로 변환
    T = table(gps_week, gps_time);

    % 테이블을 CSV 파일로 저장
    writetable(T, output_filename);
end