import numpy as np
import datetime

#输入数据
########################控制参数（由人为定义）########################
#最小沟底坡度
slope_min = 0.0008
#流速折减系数（用以计算管沟内汇流时间）
velocity_reduction_factor = 0.75
#雨力参数
rain_force = 4726
#参数n
n_factor = 0.75
#最大可接受的孔数
max_allowable_section_num = 6
#最大可接受的孔径
max_awllowable_section_width = 1.6

########################外部条件（由外部程序录入）########################
#该段沟体的所有地势高度（如路由有n段，则该列表维度为n+1）
terrain_elev = [11.68,11.69,11.54,11.44,11.54,11.44,11.35,12.28,11.35,11.59,11.53,11.44,11.54,11.49]
#沟段长度列表（如路由有n段，则该列表维度为n）
route_length = [338.5,261,162,238,162,242,115,115,34,62,146,162,128]
#上游或支线沟体接入本段排水沟的流量
subStream_flow = [0, 0, 0, 0, 0, 0, 0, 0.24, 0, 0, 0, 0, 0]
#上游或支线沟体接入本段排水沟的沟底标高（如本段路由有n段，则该列表维度为n，如汇入点在第8段沟，则在列表种第8个位置显示其沟底标高）
upstream_route_bottom_elev = [0,0,0,0,0,0,0,10.2,0,0,0,0,0]
#沟体截面形式，若为矩形，则编号为1，后期可改（如路由有n段，则该列表维度为n）
route_section_type = [1,1,1,1,1,1,1,1,1,1,1,1,1]
#沟体明暗形式，明沟为0，暗沟为1（如路由有n段，则该列表维度为n）
route_blind = [0,0,1,0,1,0,1,1,0,0,0,1,0]
#该段沟体所经区域覆土要求，如在第2至第3个点之间有一段铺筑面，铺筑面的结构总厚度为0.3，则在对应位置做输出（如路由有n段，则该列表维度为n+1）
route_coverEarth = [0,0,0.3,0.3,0,0,0,0,0,0,0,0,0,0]
#沟体顶板厚度（维度为n+1维，头一个数字为起始段的厚度）
route_coverSlab = [0.25,0.25,0.25,0.30,0.30,0.30,0.30,0.30,0.30,0.30,0.30,0.30,0.30,0.30]
#上游沟段汇入流速
initial_flow_velocity = 0
#地表汇流时间
runoff_time_on_surface = [7.05,10.51,0,6.74,0,6.33,0,0,6.4,0,8.56,0,9.83]
#累计汇水面积
accu_flowRunoff_area = [3.42,7.07,7.07,11.23,11.23,15.26,15.26,15.82,17.95,7.95,19.86,19.86,22.49]
#管沟粗糙度系数
roughness_factor = [0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013]
#沟段起始点汇入流量
route_initial_flow = 0
#初始沟段截面尺寸
section_initial_width = 1.4
section_initial_depth = 1.2
section_initial_side_slope = 0
section_num = 1
section_initial_size = [section_initial_width,section_initial_depth,section_initial_side_slope,section_num]

########################控制参数（内置参数，可由专业人员修改）########################
#当地势面坡度较小，而沟底坡度需要较大时，沟底坡度可以出现的最大值，通过观察表格而自定义的千分之五
slope_potentialHigh = 0.002
#由于规范所规定的底坡最低为0.0005，而通常排水沟穿越段顶部的最大纵坡为0.02（楼前服务车道），因此可能的解为196种(0.02-0.0005=0.0195，按步长0.0001为推算，共196种可能)，单个坡度对应的基因编码长度为6
DNA_SIZE = 6
#生成初始种群数目
POP_SIZE = 100
#代数
Iterations = 1000
#交叉率
CROSS_RATE = 0.5

########################每轮计算段落数########################
route_sec_num = len(route_length)

########################遗传算法组件（生成坡度种群并筛选）########################
class DrainageDesign_GA_calc_core():
    def __init__(self,slope_min,slope_potentialHigh,terrain_elev,route_length,route_sec_num,POP_SIZE,DNA_SIZE,Iterations,CROSS_RATE,\
                 rain_force,n_factor,route_initial_flow,max_allowable_section_num,max_awllowable_section_width,subStream_flow,\
                 route_section_type,roughness_factor,route_blind,section_initial_size,accu_flowRunoff_area,initial_flow_velocity,\
                 velocity_reduction_factor,upstream_route_bottom_elev,runoff_time_on_surface,route_coverEarth,route_coverSlab):
        #输入值
        self.slope_min = slope_min
        self.slope_potentialHigh = slope_potentialHigh
        self.terrain_elev = terrain_elev
        self.route_length = route_length
        self.route_sec_num = route_sec_num
        self.POP_SIZE = POP_SIZE
        self.DNA_SIZE = DNA_SIZE
        self.Iterations = Iterations
        self.CROSS_RATE = CROSS_RATE

        self.rain_force = rain_force
        self.n_factor = n_factor
        self.route_initial_flow = route_initial_flow
        self.max_allowable_section_num = max_allowable_section_num
        self.max_awllowable_section_width = max_awllowable_section_width
        self.subStream_flow = subStream_flow
        self.route_section_type = route_section_type
        self.roughness_factor = roughness_factor
        self.route_blind = route_blind
        self.section_initial_size = section_initial_size
        self.accu_flowRunoff_area = accu_flowRunoff_area
        self.initial_flow_velocity = initial_flow_velocity
        self.velocity_reduction_factor = velocity_reduction_factor
        self.upstream_route_bottom_elev = upstream_route_bottom_elev
        self.runoff_time_on_surface = runoff_time_on_surface
        self.route_coverEarth = route_coverEarth
        self.route_coverSlab = route_coverSlab

    #计算地形坡度，输入数组为该段沟体所投射的所有地势高度及沟段长度列表，输出为地形坡度列表(降坡为正，抬坡为负)
    def calc_terrain_slope(self):
        terrain_slope = []
        for i in range(self.route_sec_num):
            slope = round((self.terrain_elev[i] - self.terrain_elev[i+1])/self.route_length[i],4)
            terrain_slope.append(slope)
        return terrain_slope

    #确认水沟沟底坡度的上下边界
    def calc_slope_boundary(self):
        louwer_boundary = [self.slope_min]*self.route_sec_num
        terrain_slope = DrainageDesign_GA_calc_core.calc_terrain_slope(self)
        higher_boundary = []
        for i in range(self.route_sec_num):
            boundary = max(self.slope_potentialHigh, terrain_slope[i])
            higher_boundary.append(boundary)
        return [louwer_boundary,higher_boundary]

    #生成初始的坡度基因组
    def generate_DNA_group(self):
        total_dna_size = self.route_sec_num*self.DNA_SIZE
        dna_pop = np.random.randint(2, size=(self.POP_SIZE, total_dna_size))
        #预设一个拟定方案，即全沟段坡度为最小坡度
        initial_scheme = np.zeros(total_dna_size)
        #将生成的方案与拟定方案合并
        dna_pop = np.row_stack((dna_pop,initial_scheme))
        dna_pop = dna_pop[1:]
        return dna_pop

    #变异
    def crossmuta(self,dna_pop):
        new_pop = []
        for i in dna_pop:
            temp = i
            if np.random.rand() < self.CROSS_RATE:
                j = dna_pop[np.random.randint(self.POP_SIZE)]
                cpoints1 = np.random.randint(0, self.DNA_SIZE * self.route_sec_num - 1)
                cpoints2 = np.random.randint(cpoints1, DNA_SIZE * self.route_sec_num*2)
                temp[cpoints1:cpoints2] = j[cpoints1:cpoints2]
            new_pop.append(temp)
        return new_pop

    #选择
    def select(self,dna_pop,fitness):
        temp = np.random.choice(np.arange(self.POP_SIZE), size=self.POP_SIZE, replace=True, p=fitness / (fitness.sum()))
        return dna_pop[temp]

    #解码dna结果
    def decode_DNA(self,test_single_dna,route_sec_id,boundary):
        slope_scheme = test_single_dna.dot(2**np.arange(self.DNA_SIZE)[::-1])/float(2**self.DNA_SIZE-1)*(boundary[1][route_sec_id]-boundary[0][route_sec_id])+boundary[0][route_sec_id]
        return slope_scheme

    #拆解并解码DNA
    def deConstruct_DNA(self,dna_pop,boundary):
        schemes = []
        for dna_seq in dna_pop:
            slope_scheme = []
            for i in range(self.route_sec_num):
                single_dna = dna_seq[i::self.route_sec_num]
                decoded_dna = DrainageDesign_GA_calc_core.decode_DNA(self,single_dna,i,boundary)
                slope_scheme.append(round(decoded_dna,5))
            schemes.append(slope_scheme)
        return schemes

    #求解适应度
    def get_fitness(self,dna_pop):
        boundary = DrainageDesign_GA_calc_core.calc_slope_boundary(self)
        schemes = DrainageDesign_GA_calc_core.deConstruct_DNA(self, dna_pop, boundary)
        fitnessFunction = DrainageSection_fit_slope(self.rain_force, self.n_factor, self.route_initial_flow,
                                                    self.max_allowable_section_num, \
                                                    self.max_awllowable_section_width, self.subStream_flow, schemes,
                                                    self.route_sec_num, \
                                                    self.route_section_type, self.roughness_factor, self.route_blind,
                                                    self.section_initial_size, \
                                                    self.accu_flowRunoff_area, self.terrain_elev, self.route_length,
                                                    self.initial_flow_velocity, \
                                                    self.velocity_reduction_factor, self.upstream_route_bottom_elev,
                                                    self.runoff_time_on_surface, self.route_coverEarth, \
                                                    self.route_coverSlab)
        temp = fitnessFunction.calc_fitness()
        return (temp - np.min(temp)) + 0.0001

    def print_info(self,dna_pop):  # 输出最优解等
        fitness = DrainageDesign_GA_calc_core.get_fitness(self,dna_pop)
        print('fitness',fitness)
        maxfitness = np.argmax(fitness)
        print('maxfitness',maxfitness)
        print("max_fitness",fitness[maxfitness])
        print('最优的基因型:',dna_pop[maxfitness])
        boundary = DrainageDesign_GA_calc_core.calc_slope_boundary(self)
        print("解码后坡度为:",DrainageDesign_GA_calc_core.deConstruct_DNA(self,[dna_pop[maxfitness]],boundary))



    def nature_selection(self):
        start_t = datetime.datetime.now()
        dna_pop = DrainageDesign_GA_calc_core.generate_DNA_group(self)
        for i in range(Iterations):
            print('#########################第',i,'代#########################第')
            dna_pop = np.array(DrainageDesign_GA_calc_core.crossmuta(self,dna_pop))
            fitness = DrainageDesign_GA_calc_core.get_fitness(self,dna_pop)
            dna_pop = DrainageDesign_GA_calc_core.select(self,dna_pop,fitness)
        end_t = datetime.datetime.now()
        print((end_t - start_t).seconds)
        DrainageDesign_GA_calc_core.print_info(self,dna_pop)

    #输出测试数据#调试用
    def output_schemes(self):
        dna_pop = DrainageDesign_GA_calc_core.generate_DNA_group(self)
        boundary = DrainageDesign_GA_calc_core.calc_slope_boundary(self)
        schemes = DrainageDesign_GA_calc_core.deConstruct_DNA(self,dna_pop,boundary)
        return schemes

########################截面计算组件（作为遗传算法组件调用的工具，输入条件为坡度方案集合数据，输出条件为适应度值）########################
class DrainageSection_fit_slope():
    def __init__(self,rain_force,n_factor,route_initial_flow,max_allowable_section_num,max_awllowable_section_width,subStream_flow,schemes,route_sec_num,route_section_type,roughness_factor,route_blind,section_initial_size,accu_flowRunoff_area,terrain_elev,route_length,initial_flow_velocity,velocity_reduction_factor,upstream_route_bottom_elev,runoff_time_on_surface,route_coverEarth,route_coverSlab):
        self.rain_force = rain_force
        self.n_factor = n_factor
        self.route_initial_flow = route_initial_flow
        self.max_allowable_section_num = max_allowable_section_num
        self.max_awllowable_section_width = max_awllowable_section_width
        self.subStream_flow = subStream_flow
        self.schemes = schemes
        self.route_sec_num = route_sec_num
        self.route_section_type = route_section_type
        self.roughness_factor = roughness_factor
        self.route_blind = route_blind
        self.section_initial_size = section_initial_size
        self.accu_flowRunoff_area = accu_flowRunoff_area
        self.terrain_elev = terrain_elev
        self.route_length = route_length
        self.initial_flow_velocity = initial_flow_velocity
        self.velocity_reduction_factor = velocity_reduction_factor
        self.upstream_route_bottom_elev = upstream_route_bottom_elev
        self.runoff_time_on_surface = runoff_time_on_surface
        self.route_coverEarth = route_coverEarth
        self.route_coverSlab = route_coverSlab

    def calc_hydraulic_radius(self,i,section_type,section_size):
        #当沟体为矩形暗沟时计算水力半径R=(a*b)/(2*(a+b))
        if section_type == 1 and self.route_blind[i] == 1:
            hydraulic_radius = section_size[0]*section_size[1]/(2*(section_size[0]+section_size[1]))
            return hydraulic_radius
        #当沟体为梯形沟或矩形明沟时
        if section_type == 1 and self.route_blind[i] == 0:
            hydraulic_radius = (section_size[0] + section_size[1] * section_size[2]) * section_size[1] / (section_size[0] + 2 * section_size[1] * (1 + section_size[2] ** 2) ** 0.5)
            return hydraulic_radius
        #当沟体为圆管时

    def calc_section_depth(self,scheme):
        route_seg_decendent = [x * y for x, y in zip(scheme, self.route_length)]
        section_initial_depth = self.section_initial_size[1]
        initial_bottom_elev = self.terrain_elev[0]-self.route_coverEarth[0]-self.route_coverSlab[0]-section_initial_depth
        current_route_bottom_elev = [initial_bottom_elev]
        for i in range(self.route_sec_num):
            route_seg_bottom_elev = current_route_bottom_elev[-1] - route_seg_decendent[i]
            current_route_bottom_elev.append(route_seg_bottom_elev)
        current_route_cover_elev = [x-y-z for x,y,z in zip(self.terrain_elev,self.route_coverEarth,self.route_coverSlab)]
        current_route_seg_depth = [x-y for x,y in zip(current_route_cover_elev,current_route_bottom_elev)]
        return current_route_seg_depth


    def calc_section_area(self,section_type,section_size):
        # 当沟体为矩形沟时
        if section_type == 1:
            section_area = section_size[0]*section_size[1]
            return section_area

    def calc_velocity(self,i,current_slope,hydraulic_radius):
        velocity = (1/self.roughness_factor[i])*(hydraulic_radius**(2/3))*(current_slope**(1/2))
        return velocity

    def calc_runoff_time_in_ditch(self,i,velocity,current_route_upstream_velocity):
        if current_route_upstream_velocity == 0:
            runoff_time_in_current_ditch = self.route_length[i] / (60 * velocity * velocity_reduction_factor)
        else:
            runoff_time_in_current_ditch = self.route_length[i]/(30*(velocity + current_route_upstream_velocity))
        return runoff_time_in_current_ditch

    def calc_rainfall_intensity(self,rainfall_duration):
        rainfall_intensity = self.rain_force/(rainfall_duration+11.259)**self.n_factor
        return rainfall_intensity

    def compare_section_and_design_sitution(self,i,current_slope,section_size,upstream_rainfall_duration,route_seg_velocity_logger,upstream_flow):
        section_type = self.route_section_type[i]
        hydraulic_radius = DrainageSection_fit_slope.calc_hydraulic_radius(self,i,section_type, section_size)
        velocity = DrainageSection_fit_slope.calc_velocity(self,i,current_slope,hydraulic_radius)
        ###########################尚未考虑限制流速的问题##########################
        section_area = DrainageSection_fit_slope.calc_section_area(self, section_type, section_size)
        section_capacity = section_area * velocity * section_size[3]
        current_route_upstream_velocity = route_seg_velocity_logger[-1]
        runoff_time_in_current_ditch = DrainageSection_fit_slope.calc_runoff_time_in_ditch(self, i, velocity,
                                                                                           current_route_upstream_velocity)
        rainfall_duration = max(runoff_time_in_current_ditch + upstream_rainfall_duration,
                                runoff_time_in_current_ditch + self.runoff_time_on_surface[i])
        #暴雨强度
        rainfall_intensity = DrainageSection_fit_slope.calc_rainfall_intensity(self, rainfall_duration)
        route_seg_design_rainfall = max(upstream_flow,rainfall_intensity * self.accu_flowRunoff_area[i]/1000 + self.subStream_flow[i])
        diff_capacity = section_capacity - route_seg_design_rainfall
        if diff_capacity > 0:
            return [diff_capacity,runoff_time_in_current_ditch,rainfall_duration,route_seg_design_rainfall]
        else:
            return False

    def calc_fitness(self):
        #适应度记录容器
        fitness_logger = []
        #差异度记录容器
        all_scheme_diff_capacity_logger = []
        #尺寸记录容器
        all_scheme_section_logger = []
        for scheme in self.schemes:
            #分段适应度记录容器
            current_fitness_logger = []
            #分段差异度记录容器
            current_scheme_diff_capacity_logger = []
            # 截面尺寸记录容器
            scheme_size_logger = [section_initial_size]
            # 管渠内雨水流动时间
            runoff_time_in_ditch = []
            # 降雨历时记录容器
            rainfall_duration_logger = [0]
            # 截面流速记录容器
            route_seg_velocity_logger = [self.initial_flow_velocity]
            # 单跨宽度可增长次数
            width_increment_times = int((self.max_awllowable_section_width - self.section_initial_size[0]) / 0.2)
            # 设计流量记录容器
            design_flow_logger = [self.route_initial_flow]
            # 计算每段沟体的深度
            current_route_seg_depth = DrainageSection_fit_slope.calc_section_depth(self,scheme)

            for i in range(self.route_sec_num):
                #生成初始值，并确定第一段沟的初始尺寸
                section_size = scheme_size_logger[-1].copy()
                #当沟型为矩形或圆形时执行此操作
                for section_num in range(1,self.max_allowable_section_num):
                    current_slope = scheme[i]
                    section_size[1] = round(current_route_seg_depth[i+1],1)
                    section_size[3] = section_num
                    upstream_flow = design_flow_logger[-1]
                    upstream_rainfall_duration = rainfall_duration_logger[-1]
                    results = DrainageSection_fit_slope.compare_section_and_design_sitution(self,i,current_slope,section_size,upstream_rainfall_duration,route_seg_velocity_logger,upstream_flow)
                    if results:
                        scheme_size_logger.append(section_size)
                        runoff_time_in_ditch.append(results[1])
                        rainfall_duration_logger.append(results[2])
                        diff_capacity = results[0]
                        current_scheme_diff_capacity_logger.append(diff_capacity)
                        design_flow_logger.append(results[3])
                        break
                    temp_width = section_size[0]
                    for section_width_increment_times in range(1,width_increment_times):
                        section_size[0] = temp_width + section_width_increment_times * 0.2
                        section_size[0] = round(section_size[0],1)
                        results = DrainageSection_fit_slope.compare_section_and_design_sitution(self,i,current_slope,section_size,upstream_rainfall_duration,route_seg_velocity_logger,upstream_flow)
                        if results:
                            scheme_size_logger.append(section_size)
                            runoff_time_in_ditch.append(results[1])
                            rainfall_duration_logger.append(results[2])
                            diff_capacity = results[0]
                            current_scheme_diff_capacity_logger.append(diff_capacity)
                            design_flow_logger.append(results[3])
                            break
            current_scheme_diff_capacity = sum(current_scheme_diff_capacity_logger)
            current_scheme_fitness = np.log(1/current_scheme_diff_capacity)
            all_scheme_diff_capacity_logger.append(current_scheme_diff_capacity)
            all_scheme_section_logger.append(scheme_size_logger)
            fitness_logger.append(current_scheme_fitness)
        return fitness_logger


test_results = DrainageDesign_GA_calc_core(slope_min,slope_potentialHigh,terrain_elev,route_length,route_sec_num,POP_SIZE,DNA_SIZE,Iterations,CROSS_RATE,\
                 rain_force,n_factor,route_initial_flow,max_allowable_section_num,max_awllowable_section_width,subStream_flow,\
                 route_section_type,roughness_factor,route_blind,section_initial_size,accu_flowRunoff_area,initial_flow_velocity,\
                 velocity_reduction_factor,upstream_route_bottom_elev,runoff_time_on_surface,route_coverEarth,route_coverSlab)
test_results.nature_selection()































