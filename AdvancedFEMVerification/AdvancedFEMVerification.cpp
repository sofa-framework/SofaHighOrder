/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU General Public License as published by the Free  *
* Software Foundation; either version 2 of the License, or (at your option)   *
* any later version.                                                          *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    *
* more details.                                                               *
*                                                                             *
* You should have received a copy of the GNU General Public License along     *
* with this program; if not, write to the Free Software Foundation, Inc., 51  *
* Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.                   *
*******************************************************************************
*                            SOFA :: Applications                             *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/


#include <sofa/defaulttype/VecTypes.h>

//Including Simulation
#include <SofaComponentMain/init.h>
#include <sofa/simulation/graph/DAGSimulation.h>
// ForceField
#include <SofaMiscFem/StandardTetrahedralFEMForceField.h>
#include <../sofa/component/forcefield/MRForceField.h>
#include <../sofa/component/forcefield/NeoHookeanForceField.h>
#include <../sofa/component/forcefield/TetrahedralStVenantKirchhoffForceField.inl>
#include <SofaMiscFem/TetrahedralTensorMassForceField.h>
#include <SofaBaseTopology/TopologySparseData.inl>
#include <SofaBoundaryCondition/SurfacePressureForceField.h>
#include <SofaBoundaryCondition/TrianglePressureForceField.h>
#include <SofaMiscForceField/MeshMatrixMass.h>


using std::cout;
using std::cerr;
using std::endl;
using namespace sofa::component;

using namespace sofa;


const size_t sizePressureArray = 19; 

const double poissonRatioArray[] = {0.1,0.33,0.49};
const size_t sizePoissonRatioArray = sizeof(poissonRatioArray)/sizeof(poissonRatioArray[0]);

//const double pressureSVKArray[3][19]={{-.1920000000, -.1876875000, -.1785000000, -.1640625000, -.1440000000, -.1179375000, -0.8550000000e-1, -0.4631250000e-1, 0., 0.5381250000e-1, .1155000000, .1854375000, .2640000000, .3515625000, .4485000000, .5551875000, .6720000000, .7993125000, .9375000000}, {-.1920000000, -.1876875000, -.1785000000, -.1640625000, -.1440000000, -.1179375000, -0.8550000000e-1, -0.4631250000e-1, 0., 0.5381250000e-1, .1155000000, .1854375000, .2640000000, .3515625000, .4485000000, .5551875000, .6720000000, .7993125000, .9375000000}, {-.1920000000, -.1876875000, -.1785000000, -.1640625000, -.1440000000, -.1179375000, -0.8550000000e-1, -0.4631250000e-1, 0., 0.5381250000e-1, .1155000000, .1854375000, .2640000000, .3515625000, .4485000000, .5551875000, .6720000000, .7993125000, .9375000000}};
const double s1SVKArray[3][19]={{.6, .65, .70, .75, .80, .85, .90, .95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50}, {.6, .65, .70, .75, .80, .85, .90, .95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50}, {.6, .65, .70, .75, .80, .85, .90, .95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50}};
//const double s2SVKArray[3][19]={{1.031503757, 1.028469737, 1.025182911, 1.021640837, 1.017840852, 1.013780055, 1.009455299, 1.004863175, 1.000000000, .9948617996, .9894442887, .9837428525, .9777525249, .9714679617, .9648834127, .9579926931, .9507891461, .9432656044, .9354143469}, {1.100545319, 1.091134730, 1.080879272, 1.069754645, 1.057733426, 1.044784667, 1.030873416, 1.015960137, .9999999999, .9829420125, .9647279408, .9452909606, .9245539465, .9024272823, .8788060080, .8535660490, .8265591327, .7976057926, .7664854858}, {1.146123902, 1.132684864, 1.117989267, 1.101986842, 1.084619749, 1.065821280, 1.045514228, 1.023608812, 1.000000000, .9745640050, .9471536306, .9175919573, .8856635930, .8511022265, .8135723694, .7726415729, .7277362163, .6780671057, .6224949797}};


const double pressureSVKArray[3][19]={{-.1804511277, -.1774403214, -.1698382491, -.1571856287, -.1389961390, -.1147531014, -0.8390578991e-1, -0.4586531320e-1, 0., 0.5436979032e-1, .1179775280, .1916171531, .2761506277, .3725165559, .4817400644, .6049441571, .7433628317, .8983562790, 1.071428571},
{-.1585204757, -.1576444156, -.1527860993, -.1433642818, -.1287093314, -.1080434235, -0.8045544368e-1, -0.4486884492e-1, 0., 0.5569643182e-1, .1241001397, .2075231513, .3088441741, .4316960862, .5807328759, .7620183232, .9836065579, 1.256434943, 1.595744681},
{-.1461632157, -.1462908473, -.1428114248, -.1351003603, -.1224073444, -.1038205065, -0.7821791229e-1, -0.4420080647e-1, 0., 0.5665815591e-1, .1287481886, .2202410998, .3365629780, .4853321828, .6775948024, .9300012563, 1.268882174, 1.738486217, 2.419354840}};
	

const double s2SVKArray[3][19]={{1.031503757, 1.028469737, 1.025182911, 1.021640837, 1.017840852, 1.013780055, 1.009455299, 1.004863175, 1.000000000, .9948617996, .9894442887, .9837428525, .9777525249, .9714679617, .9648834127, .9579926931, .9507891461, .9432656044, .9354143469},
{1.100545319, 1.091134730, 1.080879272, 1.069754645, 1.057733426, 1.044784667, 1.030873416, 1.015960137, .9999999999, .9829420125, .9647279408, .9452909606, .9245539465, .9024272823, .8788060080, .8535660490, .8265591327, .7976057926, .7664854858},
{1.146123902, 1.132684864, 1.117989267, 1.101986842, 1.084619749, 1.065821280, 1.045514228, 1.023608812, 1.000000000, .9745640050, .9471536306, .9175919573, .8856635930, .8511022265, .8135723694, .7726415729, .7277362163, .6780671057, .6224949797}};
	


const double pressurecSVKArray[3][19]={{-.8044097165, -.6629460044, -.5474388357, -.4568349406, -.3875996268, -.3355817040, -.2968398796, -.2679192178, -.2458835407, -.2282577176, -.2129482871, -.1981688807, -.1823784835, -.1642324717, -.1425430395, -.1162436302, -0.8435129770e-1, -0.4592206473e-1, 4.545454545e-10},
{-1.409217009, -1.089788928, -.8371440806, -.6442089993, -.5004758866, -.3955378972, -.3203474723, -.2674364727, -.2307517071, -.2053928522, -.1873628056, -.1733632139, -.1606382812, -.1468588896, -.1300346984, -.1084395242, -0.8053435595e-1, -0.4487509050e-1, 0.},
{-19.77405681, -14.48953743, -10.38378450, -7.293287409, -5.019567826, -3.380185603, -2.223458495, -1.428319431, -.8995527986, -.5625393030, -.3588946719, -.2432815112, -.1812390820, -.1476699303, -.1255308429, -.1043271607, -0.7825646586e-1, -0.4420131067e-1, 1.677852333e-11}};
const double s1cSVKArray[3][19]={{.1, .15, .20, .25, .30, .35, .40, .45, .50, .55, .60, .65, .70, .75, .80, .85, .90, .95, 1.00}, {.1, .15, .20, .25, .30, .35, .40, .45, .50, .55, .60, .65, .70, .75, .80, .85, .90, .95, 1.00}, {.1, .15, .20, .25, .30, .35, .40, .45, .50, .55, .60, .65, .70, .75, .80, .85, .90, .95, 1.00}};
const double s2cSVKArray[3][19]={{1.110735356, 1.121585999, 1.123833150, 1.120397158, 1.113406234, 1.104332882, 1.094171716, 1.083588473, 1.073028292, 1.062789582, 1.053072736, 1.044011041, 1.035688721, 1.028148752, 1.021391585, 1.015365026, 1.009946880, 1.004926520, 1.000000000},
{1.204837603, 1.210401446, 1.207845760, 1.200161444, 1.189417722, 1.176976801, 1.163707695, 1.150147483, 1.136610631, 1.123259135, 1.110146037, 1.097241831, 1.084450977, 1.071624365, 1.058572609, 1.045083630, 1.030944967, 1.015966986, 1.000000000},
{1.266291029, 1.268817376, 1.263548751, 1.253513662, 1.240734771, 1.226493652, 1.211565869, 1.196388391, 1.181168424, 1.165951384, 1.150663612, 1.135142090, 1.119160501, 1.102457543, 1.084768331, 1.065853638, 1.045517891, 1.023608903, .9999999999}};

const double pressureNHArray[3][19]= { { -0.3883416901   , -0.3437710935   , -0.2998991522   , -0.2565670165   , -0.213624994   , -0.1709311133   , -0.1283498466   , -0.08575097342   , -0.04300855353   , 0.002159576726   , 0.04557678054   , 0.08950446639   , 0.1340599959   , 0.1793604752   , 0.2255234176   , 0.2726673779   , 0.3209125624   , 0.3703814535   , 0.4211994009   }, { -0.4721070502   , -0.4277836047   , -0.3835410689   , -0.3391469292   , -0.2943527318   , -0.2488897741   , -0.2024640759   , -0.1547504101   , -0.1053851729   , -0.05395775057   , 0.003699888341   , 0.08083342252   , 0.1650431622   , 0.2582399559   , 0.3628285478   , 0.4818810297   , 0.6193831983   , 0.7805911278   , 0.9725580386   }, { -0.4885595662   , -0.4468180329   , -0.4033080067   , -0.3575816025   , -0.3090534691   , -0.2569468228   , -0.200213049   , -0.1374090208   , -0.06650465361   , 0.0533210443   , 0.1127404165   , 0.1797220975   , 0.2561923322   , 0.3447263992   , 0.4488296683   , 0.5733715873   , 0.7252755222   , 0.9146526294   , 1.15673894   }};
const double s1NHArray[3][19]= { { 0.6937902788   , 0.7242478556   , 0.7555429294   , 0.7876971796   , 0.8207331994   , 0.8546745396   , 0.8895457477   , 0.9253724146   , 0.9621812234   , 1.001917966   , 1.04082847   , 1.080809799   , 1.121893624   , 1.164113005   , 1.207502465   , 1.252098068   , 1.297937492   , 1.345060128   , 1.393507161   }, { 0.5969839492   , 0.6276282091   , 0.659982157   , 0.6941793822   , 0.730367877   , 0.7687120271   , 0.8093949418   , 0.8526211953   , 0.8986200607   , 0.9476493518   , 1.003615965   , 1.079506673   , 1.162852791   , 1.254752986   , 1.356533487   , 1.469810195   , 1.596572315   , 1.739296747   , 1.901107262   }, { 0.5626500841   , 0.5926130478   , 0.6257672618   , 0.6626509616   , 0.7039306333   , 0.7504416288   , 0.8032452675   , 0.8637107801   , 0.9336357113   , 1.053287121   , 1.112386732   , 1.178304553   , 1.252292527   , 1.335929089   , 1.431233009   , 1.540828472   , 1.668190396   , 1.818019617   , 1.99683577   }};
const double s2NHArray[3][19]= { { 1.087158388   , 1.076951224   , 1.066890021   , 1.056967639   , 1.047177384   , 1.03751296   , 1.027968443   , 1.018538248   , 1.009217097   , 0.9995417907   , 0.9904288664   , 0.981410563   , 0.972482628   , 0.9636410152   , 0.9548818675   , 0.9462015031   , 0.9375964049   , 0.9290632049   , 0.9205986769   }, { 1.193240912   , 1.173969533   , 1.154716147   , 1.135469995   , 1.116220166   , 1.096955547   , 1.077664766   , 1.058336132   , 1.038957571   , 1.01951655   , 0.9986958797   , 0.9725237262   , 0.9461565385   , 0.919556827   , 0.8926836642   , 0.8654919479   , 0.837931485   , 0.8099458413   , 0.7814708799   }, { 1.323119046   , 1.290348317   , 1.256781992   , 1.222354326   , 1.186990162   , 1.150602934   , 1.113092088   , 1.074339709   , 1.034206001   , 0.9748631443   , 0.94908686   , 0.9226175014   , 0.8953935507   , 0.8673438984   , 0.8383856156   , 0.8084210025   , 0.7773336312   , 0.7449828936   , 0.71119634   }};


const double pressureMRArray[3][19]={{-.636919542677090, -.523613571530945, -.424293448456580, -.336821431383375, -.259454445400673, -.190759132592273, -.129547636959975, -0.748285075698841e-1, -0.257687595839586e-1, 0.123147654895596e-1, 0.526517764829662e-1, 0.890752787299581e-1, .122034085905884, .151915154322051, .179053559656684, .203740633528009, .226230649708071, .246746335196323, .265483452107655}, {-.562349237136935, -.493117212331550, -.427250573334016, -.364546699272354, -.304817968481825, -.247890433934346, -.193602639934889, -.141804553238043, -0.923566024350859e-1, -0.451288081083787e-1, 0.289759195812279e-1, 0.977590983272106e-1, .161678435373091, .221143906208723, .276523622983423, .328148874958982, .376318471248283, .421302493456398, .463345550710135 }, {-.507592720092027, -.456142515615975, -.404847364098504, -.353706694185020, -.302719936303602, -.251886524077834, -.201205893677865, -.150677483774722, -.100300735351388, -0.500750921217026e-1, 0.499250912708442e-1, 0.997007313853118e-1, .149327466395241, .198805839767301, .248136393195974, .297319665441725, .346356193097265, .395246510450886, .443991148046569}};
const double s1MRArray[3][19]={{.650215938710247, .687527239643059, .725929056097194, .765365159766958, .805762500759332, .847033609015900, .889079674146271, .931794051853002, .975065910151025, 1.01251520167355, 1.05652828731312, 1.10078865847721, 1.14519812021085, 1.18966661875371, 1.23411306966208, 1.27846567737549, 1.32266186026694, 1.36664786738134, 1.41037820203926}, {.639637702049833, .667909404036185, .697859697671194, .729547051724654, .763011217009673, .798268568104019, .835308062906136, .874088410079668, .914536984204839, .956550854415787, 1.02968913936779, 1.10608130132324, 1.18492093513811, 1.26540562751559, 1.34680469638069, 1.42849860927677, 1.50999193035351, 1.59090769085560, 1.67097181418148}, {.634593756547728, .659836567510541, .687211486262227, .716899878824085, .749075657156569, .783894144799765, .821478451709075, .861904401269021, .905185951795058, .951263797501961, 1.05118074082972, 1.10452763290780, 1.15971595999524, 1.21639661336373, 1.27421792985636, 1.33284425170856, 1.39196937504899, 1.45132448185498, 1.51068121696553}};
const double s2MRArray[3][19]={{1.03757670919392, 1.03394804124704, 1.02990762302508, 1.02556595581508, 1.02102439119414, 1.01637437052616, 1.01169675110348, 1.00706137054476, 1.00252693873968, .998757171436292, .994529295679769, .990512372245131, .986727246300965, .983187350362974, .979899868402788, .976866865097434, .974086327685361, .971553098458050, .969259676788529}, {1.15276938927868, 1.13801837536884, 1.12294187480164, 1.10760488960412, 1.09208196449921, 1.07645617094179, 1.06081729115063, 1.04525921482584, 1.02987669475421, 1.01476173719551, .990393133930705, .967334438704431, .945819462262607, .925976549680711, .907839716156012, .891370013246562, .876479335741178, .863051445950385, .850957995076938}, {1.24902088354625, 1.22551501674744, 1.20146363702156, 1.17691573903082, 1.15194130295011, 1.12663393953631, 1.10111221710127, 1.07551885873374, 1.05001709087481, 1.02478381944566, .975839382171100, .952457416859653, .929982228608163, .908509044567281, .888098578545382, .868778948075975, .850550090414534, .833389463544910, .817257977203409}};


const double pressureMR2Array[3][19]={{-0.375000000294828, -0.331250000363882, -0.287500000064940, -0.243750000188956, -0.200000000553883, -0.156250000021032, -0.112500000268643, -0.687500004759214e-1, -0.250000003727931e-1, 0.125000004280802e-1, 0.562499992531911e-1, 0.999999993289800e-1,0.143750000192142,0.187499999718972,0.231249999353961,0.274999999858699,0.318749999527915,0.362500000107645,0.406250000602761}, {-0.441176471538028, -0.397058824685343, -0.352941175994425, -0.308823529665840, -0.264705882542141, -0.220588234457488, -0.176470588496682, -0.132352940737952, -0.882352944670396e-1, -0.441176462788487e-1, 0.294117658219286e-1,0.102941177746872,0.176470588078132,0.250000000690293,0.323529413221062,0.397058824022121,0.470588235448882,0.544117646626235,0.617647058086856}, {-0.500000016167740, -0.450000027576502, -0.400000010529129, -0.350000010086341, -0.300000011685702, -0.250000000157779, -0.199999991216094, -0.150000000946984, -0.999999847869861e-1, -0.499999891465964e-1, 0.499999951887867e-1,0.100000010717559,0.149999999058132,0.200000003072519,0.250000005365214,0.300000047767465,0.349999994699014,0.399999976046932,0.450000002884075}};
const double s1MR2Array[3][19]={{0.734923710696966,0.756964104906148,0.780832188656682,0.806745032705690,0.834935315042303,0.865650942857756,0.899153173176482,0.935712861720541,0.975604496086924, 1.01265410800356, 1.05943801535170, 1.11027453412948, 1.16535629510163, 1.22482486190523, 1.28876098048364, 1.35717896744314, 1.43002632501880, 1.50718874305059, 1.58849962930757}, {.690356977657755,0.711598578715659,0.734562547900050,0.759406973239709,0.786298284315693,0.815408389734820,0.846910514084144,0.880973608981182,0.917755352628016,0.957394052027525, 1.03009070915148, 1.11128788359967, 1.20090971526327, 1.29852330484387, 1.40340405916439, 1.51465123000230, 1.63131260317050, 1.75248381813230, 1.87736799637715}, {.638191049142679,0.662995883673803,0.689919182697698,0.719148058442466,0.750865116377410,0.785237909899124,0.822405719629479,0.862464534238259,0.905451964210824,0.951334434830444, 1.05125874630161, 1.10485192442050, 1.16046855886543, 1.21776691875409, 1.27639651420617, 1.33601782872502, 1.39631688138978, 1.45701509831430, 1.51787314976536}};
const double s2MR2Array[3][19]={{0.975950723650987,0.985385180782041,0.993039465208973,0.998917887742305, 1.00302839485190, 1.00538551369287, 1.00601360848852, 1.00495034755551, 1.00225017625436,0.998688669047359,0.993162627966488,0.986271965014372,0.978155881093731,0.968973552792587,0.958899265320095,0.948115736729822,0.936806394981441,0.925147616471180,0.913301940143439}, {1.10961572686937, 1.10253025950479, 1.09452817411076, 1.08561098029022, 1.07578883366999, 1.06508257251996, 1.05352577319253, 1.04116662094902, 1.02806933206304, 1.01431477747675,0.990200068057227,0.965065710901562,0.939502096376426,0.914092182726047,0.889344708131952,0.865648557831943,0.843257861835180,0.822303553291127,0.802819624578315}, {1.24549572987470, 1.22259161028792, 1.19910365262723, 1.17507467926404, 1.15056783282953, 1.12566953142825, 1.10049128635074, 1.07516955079513, 1.04986283746759, 1.02474577323794,0.975803176896283,0.952317625943857,0.929680618933709,0.907997745878644,0.887340340827748,0.867746485740476,0.849224941821046,0.831760400288946,0.815319526901478}};

const double pressureBAArray[3][19]= { { -0.828428684899   , -0.672216813836   , -0.538888261552   , -0.424293448219   , -0.325186776574   , -0.238998425601   , -0.163671090341   , -0.0975412536752   , -0.039251555945   , 0.0012481273851   , 0.0697594954508   , 0.127365236562   , 0.176081220653   , 0.217483386676   , 0.252818216369   , 0.283082038521   , 0.309078861034   , 0.331463145978   , 0.350771866002   }, { -0.711803940807   , -0.635165741897   , -0.562349236902   , -0.493117212129   , -0.427250572965   , -0.364546699722   , -0.304817968625   , -0.24789043392   , -0.193602639858   , -0.141804553435   , -0.0923566023925   , -0.045128808038   , 0.00293677027376   , 0.059867889446   , 0.113527122273   , 0.164139988714   , 0.211913094076   , 0.257036003874   , 0.299682909504   }, { -0.662879410219   , -0.57471073524   , -0.486993994815   , -0.399726354642   , -0.312905001904   , -0.226527146129   , -0.140590018743   , -0.0550908721949   , 0.0499250915536   , 0.0997007317344   , 0.149327466563   , 0.198805840212   , 0.248136393484   , 0.297319665859   , 0.346356193559   , 0.395246510746   , 0.443991148374   , 0.492590636229   , 0.541045501941   }};
const double s1BAArray[3][19]= { { 0.57243167835   , 0.620561394356   , 0.669568582231   , 0.719087658315   , 0.768786090503   , 0.81837821737   , 0.867630697388   , 0.916361547277   , 0.964435169115   , 1.00117782895   , 1.07089808729   , 1.13868181379   , 1.2044904493   , 1.26834452957   , 1.330300332   , 1.3904342203   , 1.44883247319   , 1.50558484992   , 1.56078063972   }, { 0.540881270878   , 0.572696513167   , 0.606240159453   , 0.641407484823   , 0.678052141816   , 0.71599241638   , 0.755020852734   , 0.794915862744   , 0.835453607907   , 0.876418587941   , 0.917611859953   , 0.958856440149   , 1.00273605544   , 1.0571936569   , 1.11100198248   , 1.1639937029   , 1.21605533453   , 1.2671158196   , 1.31713642615   }, { 0.501787443497   , 0.54779782441   , 0.600409565393   , 0.659781378525   , 0.725535442228   , 0.796700665777   , 0.871841617879   , 0.94933347437   , 1.04605683367   , 1.09188010847   , 1.13729020442   , 1.1821466373   , 1.2263434657   , 1.26980403445   , 1.31247575787   , 1.35432534862   , 1.39533467649   , 1.4354973028   , 1.47481565051   }};
const double s2BAArray[3][19]= { { 1.06560156222   , 1.05446499783   , 1.04415241286   , 1.03479529421   , 1.02645461312   , 1.01913667889   , 1.0128089925   , 1.00741414824   , 1.00288087101   , 0.999911186186   , 0.995367876174   , 0.992205463388   , 0.990200789917   , 0.989161650826   , 0.988926007753   , 0.989358840322   , 0.990348283669   , 0.991801843672   , 0.993643016743   }, { 1.23127753156   , 1.20748284203   , 1.18409653917   , 1.16129095101   , 1.13922574268   , 1.11803990971   , 1.09784549205   , 1.07872385932   , 1.06072488424   , 1.04386874855   , 1.02814971172   , 1.01354099747   , 0.999133965764   , 0.982733128645   , 0.96798348508   , 0.954733972734   , 0.942836047052   , 0.932149354108   , 0.922544716621   }, { 1.40248647429   , 1.34345188035   , 1.28434445355   , 1.22624626569   , 1.17036191256   , 1.11782300218   , 1.06947918561   , 1.02577381636   , 0.978226442879   , 0.957957812608   , 0.93910643669   , 0.921576052272   , 0.905267585756   , 0.890083271944   , 0.875929373029   , 0.862717807523   , 0.850366982721   , 0.838802075659   , 0.827954953868   }};


/**  Test force fields implementing linear elasticity on tetrahedral mesh. 
Implement traction applied on the top part of a cylinder and test that the deformation 
is simply related with the Young Modulus and Poisson Ratio of the isotropc linear elastic material */


    typedef  sofa::defaulttype::Vec3dTypes DataTypes;
    typedef  DataTypes::VecCoord VecCoord;
    typedef  DataTypes::VecDeriv VecDeriv;
    typedef  DataTypes::Deriv Deriv;
	typedef  sofa::defaulttype::Vec3d Coord;
	typedef  DataTypes::Real Real;
	typedef const double dataArray[3][19];
    typedef  container::MechanicalObject<DataTypes> MechanicalObject;
    typedef  sofa::component::mass::MeshMatrixMass<DataTypes,Real>  MeshMatrixMass;
	typedef  sofa::component::forcefield::StandardTetrahedralFEMForceField<DataTypes> StandardTetrahedralFEMForceField;
	typedef  sofa::component::forcefield::MRForceField<DataTypes> MooneyRivlinForceField;
	typedef  sofa::component::forcefield::NeoHookeanForceField<DataTypes> NeoHookeanAFForceField;
	typedef  sofa::component::forcefield::TetrahedralStVenantKirchhoffForceField<DataTypes> TetrahedralSVKForceField;
    typedef  sofa::core::behavior::ForceField<DataTypes>::SPtr ForceFieldSPtr;
    typedef ForceFieldSPtr (*HyperElasticFF)(simulation::Node::SPtr,double,double);
	typedef  sofa::component::forcefield::TrianglePressureForceField<DataTypes> TrianglePressureForceField;
	typedef  sofa::component::forcefield::SurfacePressureForceField<DataTypes> SurfacePressureForceField;
	
    /// Simulation
    simulation::Simulation* simulation;
	simulation::Node::SPtr root;
	/// struct with the pointer of the main components 
//	CylinderTractionStruct<DataTypes> tractionStruct;
	/// index of the vertex used to compute the compute the deformation
	size_t vIndex;
	sofa::helper::vector<size_t> indices;


     // Define the path for the scenes directory
    #define ADD_SOFA_TEST_SCENES_PATH( x ) sofa_tostring(SOFA_ADVANCEDFEM_VERIFICATION_SCENES_PATH)sofa_tostring(x) 
	

/** Helper class to create a component and add it as a child of a given Node */
template<class T>
class addNew : public core::objectmodel::New<T>
{
    typedef typename T::SPtr SPtr;
public:
    addNew( simulation::Node::SPtr parent, const char* name="")
    {
        parent->addObject(*this);
        (*this)->setName(name);
    }
};


	 void loadScene(std::string sceneName)
    {
        // Get the scene directory
        sofa::helper::system::FileRepository repository("SOFA_DATA_PATH");
        repository.addFirstPath( ADD_SOFA_TEST_SCENES_PATH( /Scenes ) );


        // Load the scene from the xml file
        std::string fileName = repository.getFile(sceneName);
		/// Root of the scene graph
	   root = sofa::core::objectmodel::SPtr_dynamic_cast<sofa::simulation::Node>( sofa::simulation::getSimulation()->load(fileName.c_str()));
		
	 }
    // Create the context for the scene
    void SetUp()
    { 
        // Init simulation
        sofa::component::init();
		sofa::simulation::graph::DAGSimulation *simu = new sofa::simulation::graph::DAGSimulation();
        sofa::simulation::setSimulation(simu);
	}
   void LoadTetraTraction()
   {
		loadScene("HyperelasticTestScene.scn");	 

		size_t resolutionCircumferential=7;
		size_t  resolutionRadial=3;
		size_t  resolutionHeight=7;
	
		/// take the vertex at mid height and on the surface of the cylinder
		vIndex=(resolutionCircumferential*(resolutionRadial-1)+1)*resolutionHeight/2;
    }
      void LoadTetraSimpleShear()
   {
		loadScene("HyperelasticSimpleShearTetraTestScene.xml");	 
		/// add the indices of the 8 corners of the cube
		indices.push_back(0);
		indices.push_back(5);
		indices.push_back(30);
		indices.push_back(35);
		indices.push_back(180);
		indices.push_back(185);
		indices.push_back(210);
		indices.push_back(215);
	

    }
	ForceFieldSPtr addTetrahedralSVKMaterial(simulation::Node::SPtr root,
        double youngModulus,double poissonRatio)
	{
		TetrahedralSVKForceField::SPtr ff=addNew<TetrahedralSVKForceField>(root);
		ff->setYoungModulus(youngModulus);
		ff->setPoissonRatio(poissonRatio);

		std::string str;
		ff->f_compressible.setValue(false);
		
		return (ForceFieldSPtr )ff; 
	}
	ForceFieldSPtr addTetrahedralStandardSVKMaterial(simulation::Node::SPtr root,
        double youngModulus,double poissonRatio)
	{
        StandardTetrahedralFEMForceField::SPtr ff=addNew<StandardTetrahedralFEMForceField>(root);
		Real lambda=youngModulus*poissonRatio/((1+poissonRatio)*(1-2*poissonRatio));
		Real mu=youngModulus/(2*(1+poissonRatio));
		std::vector<Real> pSet(2);
		pSet[1]=lambda; pSet[0]=mu;
		ff->setparameter(pSet);
        ff->setMaterialName("StVenantKirchhoff");
		return (ForceFieldSPtr )ff;
	}
	ForceFieldSPtr addTetrahedralcSVKMaterial(simulation::Node::SPtr root,
        double youngModulus,double poissonRatio)
	{
         TetrahedralSVKForceField::SPtr ff=addNew<TetrahedralSVKForceField>(root);
		ff->setYoungModulus(youngModulus);
		ff->setPoissonRatio(poissonRatio);
		ff->f_compressible.setValue(true);
		Real lambda=youngModulus*poissonRatio/((1+poissonRatio)*(1-2*poissonRatio));
		Real mu=youngModulus/(2*(1+poissonRatio));
		return (ForceFieldSPtr )ff; 
	}
	ForceFieldSPtr addTetrahedralStandardNeoHookeanMaterial(simulation::Node::SPtr root,
        double youngModulus,double poissonRatio)
	{
        StandardTetrahedralFEMForceField::SPtr ff=addNew<StandardTetrahedralFEMForceField>(root);
		Real lambda=youngModulus*poissonRatio/((1+poissonRatio)*(1-2*poissonRatio));
		Real mu=youngModulus/(2*(1+poissonRatio));
		Real bulkModulus=lambda+2*mu/3;
		sofa::helper::vector<Real> pSet(2);
		pSet[1]=bulkModulus; pSet[0]=mu;
		std::cerr<<"material= "<<pSet<<std::endl;
		ff->setparameter(pSet);
        ff->setMaterialName("NeoHookean");
		return (ForceFieldSPtr )ff;
	}
ForceFieldSPtr addTetrahedralStandardBoyceArrudaMaterial(simulation::Node::SPtr root,
        double youngModulus,double poissonRatio)
	{
        StandardTetrahedralFEMForceField::SPtr ff=addNew<StandardTetrahedralFEMForceField>(root);
		Real lambda=youngModulus*poissonRatio/((1+poissonRatio)*(1-2*poissonRatio));
		Real mu=youngModulus/(2*(1+poissonRatio));
		Real bulkModulus=lambda+2*mu/3;
		sofa::helper::vector<Real> pSet(2);
		pSet[1]=bulkModulus; pSet[0]=mu;
		ff->setparameter(pSet);
        ff->setMaterialName("ArrudaBoyce");
		return (ForceFieldSPtr )ff;
	}
ForceFieldSPtr addTetrahedralAdvancedFEMNeoHookeanMaterial(simulation::Node::SPtr root,
        double youngModulus,double poissonRatio)
	{
        NeoHookeanAFForceField::SPtr ff=addNew<NeoHookeanAFForceField>(root);
		Real lambda=youngModulus*poissonRatio/((1+poissonRatio)*(1-2*poissonRatio));
		Real mu=youngModulus/(2*(1+poissonRatio));
		Real bulkModulus=lambda+2*mu/3;
		ff->f_parametermu.setValue(mu);
		ff->f_parameterko.setValue(bulkModulus);
		return (ForceFieldSPtr )ff;
	}
ForceFieldSPtr addTetrahedralStandardMooneyRivlinMaterial(simulation::Node::SPtr root,
        double youngModulus,double poissonRatio)
	{
        StandardTetrahedralFEMForceField::SPtr ff=addNew<StandardTetrahedralFEMForceField>(root);
		Real lambda=youngModulus*poissonRatio/((1+poissonRatio)*(1-2*poissonRatio));
		Real mu=youngModulus/(2*(1+poissonRatio));
		Real bulkModulus=lambda+2*mu/3;
		std::vector<Real> pSet(3);
		pSet[2]=bulkModulus; pSet[0]=mu/4; pSet[1]=mu/4;
		ff->setparameter(pSet);
        ff->setMaterialName("MooneyRivlin");
		return (ForceFieldSPtr )ff;
	}
ForceFieldSPtr addTetrahedralAdvancedFEMMooneyRivlinMaterial(simulation::Node::SPtr root,
        double youngModulus,double poissonRatio)
	{
        MooneyRivlinForceField::SPtr ff=addNew<MooneyRivlinForceField>(root);
//		 typename MooneyRivlinForceField *ff = new MooneyRivlinForceField
		Real lambda=youngModulus*poissonRatio/((1+poissonRatio)*(1-2*poissonRatio));
		Real mu=youngModulus/(2*(1+poissonRatio));
		Real bulkModulus=lambda+2*mu/3;
		Coord pSet;
		pSet[2]=bulkModulus; pSet[0]=mu/4; pSet[1]=mu/4;
		helper::vector<Coord > *pps=ff->f_parameterSet.beginEdit();
		(*pps).push_back(pSet);
		ff->f_parameterSet.endEdit();
//		ff->f_parameterSet.push_back(pSet);
		return (ForceFieldSPtr )ff;
	}

	bool testHyperElasticMaterialInSimpleShearWithTetras(HyperElasticFF createForceField, dataArray pressureArray, 
		dataArray s1Array,dataArray s2Array,dataArray s3Array,size_t minIdex,size_t maxIndex){
					
		size_t i,k;
		std::vector<void *> listTPF;
//		TrianglePressureForceField *tpf;
		root->getContext()->get<TrianglePressureForceField,std::vector<void *> >(&listTPF);
		if (listTPF.size()!=4) {
			cerr << "Could not find the 4 Pressure Force fields "<< std::endl;
			exit(1);
		}
		MechanicalObject *dof;
		root->getContext()->get(dof);
		if (!dof) {
			std::cerr << "Could not find Mechanical Object "<< std::endl;
			exit(1);
		}
		Real youngModulus=1.0;
		Real poissonRatio=poissonRatioArray[0];

		
			ForceFieldSPtr ff=(*createForceField)(root,youngModulus,0.45);
			ff->init();
			// set the shear stress of the 4 faces of the cube 
			Real pressure= pressureArray[0][0];
			pressure= 0.1244961005;
			for (i=0;i<4;++i) {
				static_cast<TrianglePressureForceField *>(listTPF[i])->cauchyStress.setValue(defaulttype::MatSym<3,Real>(0.0f,0.0f,0.0f,pressure,0.0f,0.0f));
				static_cast<TrianglePressureForceField *>(listTPF[i])->p_useConstantForce.setValue(false);
			}
			TrianglePressureForceField *tpf0=static_cast<TrianglePressureForceField *>(listTPF[0]);
			tpf0->init();
			// reset simulation and init the triangle pressure forcefield
				sofa::simulation::getSimulation()->reset(root.get());
				for (i=0;i<4;++i)
					static_cast<TrianglePressureForceField *>(listTPF[i])->init();
				
				// record the initial point of 8 vertices
				Coord p[8],q[8];
				for (i=0;i<8;++i)
					p[i]= dof->read(core::ConstVecCoordId::position())->getValue()[indices[i]];
				
				//  do several steps of the static solver
				for(k=0;k<20;++k) {
					sofa::simulation::getSimulation()->animate(root.get());
				}
				// record the final position of 8 vertices
				for (i=0;i<8;++i) {
					q[i]= dof->read(core::ConstVecCoordId::position())->getValue()[indices[i]];
						std::cerr<< "New pos ="<<q[i]<<" old pos="<<p[i]<<std::endl;
				}
				// check that it is a parallepiped
				if (((q[7]-q[6]-q[1]+q[0]).norm()>1e-4) ||  ((q[7]-q[5]-q[2]+q[0]).norm()>1e-4) ||  ((q[7]-q[3]-q[4]+q[0]).norm()>1e-4)){
					std::cerr<< "Not an affine transformation q[7]-q[6]-q[1]+q[0]).norm()="<<(q[7]-q[6]-q[1]+q[0]).norm()   <<" (q[7]-q[5]-q[2]+q[0]).norm()old pos="
						<<(q[7]-q[5]-q[2]+q[0]).norm()<< " (q[7]-q[3]-q[4]+q[0]).norm()="<<(q[7]-q[3]-q[4]+q[0]).norm() <<std::endl;
					return (false);
				}
				Real s1=q[1][0];
				Real s2=q[2][1];
				Real s3=q[4][2];
				std::cerr<< " s1= "<<s1<<" s2="<<s2<<" s3="<<s3<<std::endl;
				if ((fabs(q[7][0]-s1-sqrt(s3*s3-s1*s1))>1e-4) || (fabs(q[7][1]-s2)>1e-4)|| (fabs(q[7][2]-s3)>1e-4)){
					std::cerr<< "Error in deformation since  fabs(q[7][0]-s1-sqrt(1-s1*s1/s2*s2)*s3)="<<fabs(q[7][0]-s1-sqrt(1-s1*s1/s2*s2)*s3)<< " and fabs(q[7][1]-s2)="<<
						fabs(q[7][1]-s2)<< " and fabs(q[7][2]-s3)="<<fabs(q[7][2]-s3)<<std::endl;
					return (false);
				}
				std::cerr<< " s1= "<<s1<<" s2="<<s2<<" s3="<<s3<<std::endl;
				return true;
	}
	bool testHyperElasticMaterialInTraction(HyperElasticFF createForceField, dataArray pressureArray, 
		dataArray s1Array,dataArray s2Array,size_t minIdex,size_t maxIndex){
					
		
		

		SurfacePressureForceField *tpf;
		root->getContext()->get(tpf);
		if (!tpf) {
			cerr << "Could not find Pressure Force field "<< std::endl;
			exit(1);
		}
		MechanicalObject *dof;
		root->getContext()->get(dof);
		if (!dof) {
			std::cerr << "Could not find Mechanical Object "<< std::endl;
			exit(1);
		}
		size_t i,j,k;
		Real youngModulus=1.0;
		for (j=0;j<sizePoissonRatioArray;++j) {
			Real poissonRatio=poissonRatioArray[j];
			// create the linear elasticity force field
			
			ForceFieldSPtr ff=(*createForceField)(root,youngModulus,poissonRatio);
			ff->init();
	
			for (i=minIdex;i<maxIndex;++i) {
				// set the pressure on the top part
				Real pressure= pressureArray[j][i];
				//tpf->pressure= Coord(0,0,pressure);
				tpf->setPressure(-pressure);
	//			tpf->pressureScalar= -pressure;

				// reset simulation and init the triangle pressure forcefield
				sofa::simulation::getSimulation()->reset(root.get());
				// sofa::simulation::getSimulation()->init(tractionStruct.root.get());
				tpf->init();
				// record the initial point of a given vertex
                Coord p0 = dof->read(core::ConstVecCoordId::position())->getValue()[vIndex];

				//  do several steps of the static solver
				for(k=0;k<20;++k) {
					sofa::simulation::getSimulation()->animate(root.get(),0.5);
				}
	//			std::cerr << "pressure no "<<i<< "with value "<<  pressure <<std::endl;
				// Get the simulated final position of that vertex
                Coord p1 = dof->read(core::ConstVecCoordId::position())->getValue()[vIndex];
				if (p1[0]!=p1[0]) {
					std::cerr << "Evaluation not converging for Poisson Ratio = "<<
						poissonRatio << " pressure= "<<pressure<< std::endl;
					return false;
				}
				// test the longitudinal stretch
				Real longitudinalStretch=p1[2]/p0[2];
				if (fabs(longitudinalStretch-s1Array[j][i])>1e-3) {
					std::cerr << "Wrong longitudinal stretch for Poisson Ratio = "<<
						poissonRatio << " pressure= "<<pressure<< std::endl <<
						"Got "<<longitudinalStretch<< " instead of "<< s1Array[j][i]<< std::endl;
					return false;
				}
				// compute radial stretch
				p0[2]=0;
				p1[2]=0;
				Real radius=p0.norm2();
				Real radialStretch= dot(p0,p1)/radius;
				// test the radial stretch
				if (fabs(radialStretch-s2Array[j][i])>1e-3) {
					std::cerr<< "Wrong radial stretch for  Poisson Ratio = "<<
						poissonRatio << " pressure= "<<pressure<< std::endl <<
						"Got "<<radialStretch<< " instead of "<< s2Array[j][i]<< std::endl;
					return false;
				}
			}
			root->removeObject(ff);
		}

		return true;
	}
void testTetrasInSimpleShear() {
	LoadTetraSimpleShear();
	sofa::simulation::getSimulation()->init(root.get());
	if (testHyperElasticMaterialInSimpleShearWithTetras(&addTetrahedralStandardNeoHookeanMaterial,
		pressureNHArray,s1NHArray, s2NHArray,s2NHArray,7,19 ))
	{
		std::cerr<<"Passed Tests on sheared TetrahedralFEMForceField"<<std::endl;
	}
	exit(1);
}

void testTetrasInTraction() {
	LoadTetraTraction();
	sofa::simulation::getSimulation()->init(root.get());

	if (testHyperElasticMaterialInTraction(&addTetrahedralcSVKMaterial,
		pressurecSVKArray,s1cSVKArray, s2cSVKArray,7,19 ))
	{
		std::cerr<<"Passed Tests on compressed TetrahedralFEMForceField"<<std::endl;
	}
	
	sofa::simulation::getSimulation()->reset(root.get());
	if (testHyperElasticMaterialInTraction(&addTetrahedralStandardSVKMaterial,
		pressureSVKArray,s1SVKArray, s2SVKArray,3,19 ))
	{
		std::cerr<<"Passed Tests on standard SVK"<<std::endl;
	}
	sofa::simulation::getSimulation()->reset(root.get());

	if (testHyperElasticMaterialInTraction(&addTetrahedralSVKMaterial,
		pressureSVKArray,s1SVKArray, s2SVKArray,2,19 ))
	{
		std::cerr<<"Passed Tests on TetrahedralSVKForceField"<<std::endl;
	} 
	sofa::simulation::getSimulation()->reset(root.get());

	if (testHyperElasticMaterialInTraction(&addTetrahedralStandardMooneyRivlinMaterial,
		pressureMRArray,s1MRArray, s2MRArray,0,19 ))
	{
		std::cerr<<"Passed Tests on standard Moonley Rivlin"<<std::endl;
	} 
	sofa::simulation::getSimulation()->reset(root.get());
	if (testHyperElasticMaterialInTraction(&addTetrahedralAdvancedFEMMooneyRivlinMaterial,
		pressureMRArray,s1MRArray, s2MRArray,2,19 ))
	{
		std::cerr<<"Passed Tests on Advanced Moonley Rivlin"<<std::endl;
	} 
	sofa::simulation::getSimulation()->reset(root.get());

	if (testHyperElasticMaterialInTraction(&addTetrahedralStandardNeoHookeanMaterial,
		pressureNHArray,s1NHArray, s2NHArray,0,19 ))
	{
		std::cerr<<"Passed Tests on standard Neo Hookean"<<std::endl;
	} 
	sofa::simulation::getSimulation()->reset(root.get());

	if (testHyperElasticMaterialInTraction(&addTetrahedralAdvancedFEMNeoHookeanMaterial,
		pressureNHArray,s1NHArray, s2NHArray,0,19 ))
	{
		std::cerr<<"Passed Tests on advanced  Neo Hookean"<<std::endl;
	} 
	sofa::simulation::getSimulation()->reset(root.get());

	if (testHyperElasticMaterialInTraction(&addTetrahedralStandardBoyceArrudaMaterial,
		pressureBAArray,s1BAArray, s2BAArray,1,19 ))
	{
		std::cerr<<"Passed Tests on standard Arruda Boyce"<<std::endl;
	} 
	sofa::simulation::getSimulation()->reset(root.get());	
	
        sofa::simulation::getSimulation()->unload(root.get());
}

int main(int argc, char** argv)
{
	SetUp();
	testTetrasInSimpleShear();
	testTetrasInTraction();

	
}



