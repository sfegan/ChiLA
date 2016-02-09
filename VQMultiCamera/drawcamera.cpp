//-*-mode:c++; mode:font-lock;-*-

#include <fstream>
#include <sstream>
#include <string>
#include <cctype>
#include <cmath>

#include <Qt/QtGui>

#include <RandomNumbers.hpp>
#include <VSOptions.hpp>
#include <VSLineTokenizer.hpp>

#include "VSimpleCameraLayout.hpp"
#include "VQMultiCamera.hpp"
#include "VQColormap.hpp"

using namespace VERITAS;

static std::string toLower(const std::string& s)
{
  std::string t;
  for(unsigned i=0; i<s.length(); i++)t += tolower(s[i]);
  return t;
}

static std::string toUpper(const std::string& s)
{
  std::string t;
  for(unsigned i=0; i<s.length(); i++)t += toupper(s[i]);
  return t;
}

void 
blendColors(qreal x, QColor& c, const QColor& b)
{
  if(x<0)x=0;
  if(x>1)x=1;
  qreal omx=1-x;
  QColor temp;
  temp.setRgbF(x*c.redF()   + omx*b.redF(),
	       x*c.greenF() + omx*b.greenF(),
	       x*c.blueF()  + omx*b.blueF());
  c=temp;
}

void usage(const std::string& progname, const VSOptions& options, 
	   std::ostream& stream)
{
  stream << "Usage: " << progname << " [options] [filename]" << std::endl
	 << std::endl;
  stream << "Options:" << std::endl;
  options.printUsage(stream);
}

int main(int argc, char** argv)
{
  bool batch = false;
  for(int iarg=1;iarg<argc;iarg++)
    {
      if((strcmp(argv[iarg],"-batch")==0)
	 ||(strcmp(argv[iarg],"--batch")==0)
	 ||(strncmp(argv[iarg],"-batch=",7)==0)
	 ||(strncmp(argv[iarg],"--batch=",8)==0))
	batch=true;
    }

  QApplication app(argc,argv,!batch);
  int my_argc = app.argc();
  char** my_argv = new char*[my_argc+1];
  for(int iarg=0;iarg<my_argc;iarg++)my_argv[iarg]=app.argv()[iarg];
  my_argv[my_argc]=NULL;

  VSOptions options(my_argc, my_argv, true);

  bool draw_fast = false;
  std::string camera_type = "STANDARD";
  double tube_diameter = 0.074*2;
  double camera_diameter = 1.666*2;
  unsigned hex_rings = 12;
  std::string whipple_camera = "VC500BIGTUBES";
  std::string draw_type = "CIRCLE";
  bool no_draw_outline = false;
  std::string png_file = "";
  std::string ps_file = "";
#if ( QT_VERSION >= 0x040100 )
  std::string pdf_file = "";
#endif
  std::string value_mode = "color";
  std::string color_mode = "scope";
  std::string colormap = "blue_red_yellow";
  std::vector<std::string> colors;
  std::vector<std::string> columns;
  bool log_scale = false;
  bool subtract_mean = false;
  std::pair<double,double> limits = std::make_pair<double,double>(0,0);
  bool clip_data = false;
  double value_default = 0;
  unsigned value_column = 0;
  unsigned png_size = 2048;
  unsigned ncam = 1;
  bool camera_overlap = false;

  bool label_channels = false;
  std::string font_family = "helvetica";
  double font_scale = 0.9;

  bool print_usage = false;
  bool no_pause = false;

  if(options.find("h","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;
  if(options.find("help","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;
  
  options.find("batch","Batch mode: disable GUI.");

  if(options.find("no_pause","Do not pause waiting for the QT window to be "
		  "closed. This option should not be necessary since the "
		  "\"batch\" option is cleaner. However there seems to be "
		  "some problems with \"batch\" under some QT versions.")
     !=VSOptions::FS_NOT_FOUND)
    no_pause=true;

  if(options.find("fast",
		  "Draw X11 image in \"fast\" mode, sacrificing accuracy")
     !=VSOptions::FS_NOT_FOUND)
    draw_fast=true;

  options.findWithValue("camera_type", camera_type,
			"The type of camera to draw. Valid options are: "
			"(1) \"STANDARD\", draw one of the standard "
			"Whipple/VERITAS cameras specified by the "
			"\"-camera\" option, (2) \"ROUNDED\" draw a rounded "
			"camera of a given diameter or (3) \"HEX\" draw a "
			"hexagonal radius with a given number of rings.");

  options.findWithValue("camera", whipple_camera,
			"The type of standard camera to draw, if "
			"\"-camera_type=STANDARD\". Valid options are: "
			"WC109, WC151, WC331, WC490, VC499, VC499BIGTUBES, "
			"VC500 and VC500BIGTUBES.");

  options.findWithValue("camera_tube_size", tube_diameter,
			"Set the tube diameter for the \"ROUNDED\" and "
			"\"HEX\" camera types.");

  options.findWithValue("camera_diameter", camera_diameter,
			"Set the camera diameter for the \"ROUNDED\" "
			"camera type.");

  options.findWithValue("camera_rings", hex_rings,
			"Set the number of hexagonal rings in the \"HEX\" "
			"camera type. The number of channels is 3*N*(N+1)+1.");

  options.findWithValue("camera_num", ncam,
			"Set the number of cameras to draw");

  if(options.find("camera_overlap",
		  "Do not draw multiple seperated cameras, instead draw all "
		  "camera data overlapping")!=VSOptions::FS_NOT_FOUND)
    camera_overlap=true;

  options.findWithValue("tube_type", draw_type,
			"Set the type (shape) of tube to draw. Choices are: "
			"\"CIRCLE\", \"SQUARE\", \"SMALLSQUARE\", \"HEX\", "
			"\"SMALLHEX\", \"STAR\", \"PIE\".");

  options.findWithValue("tube_value", value_mode, 
			"Set the mode to use for expressing the data values. "
			"Valid options are: \"area\", draw the area of the "
			"channel in proportion to the value; \"radius\", draw "
			"the radius in proportion to the value and \"color\", "
			"use a color map to color the channel in proportion "
			"to the value.");

  options.findWithValue("tube_color", color_mode,
			"Set the method of coloring the channel, if "
			"\"value_mode\" is NOT \"color\". Valid options are: "
			"\"scope\", channels are colored by telescope ID or "
			"\"type\", channels are colored by data type\".");

  if(options.find("tube_label", "Enable labeling of channels.") 
     != VSOptions::FS_NOT_FOUND)
    label_channels = true;

  if(options.find("tube_disable_outline", "Disable drawing of tube outline")
     !=VSOptions::FS_NOT_FOUND)
    no_draw_outline=true;
  
  options.findWithValue("font_family",font_family,
			"Set the font family for labeling channels.");

  options.findWithValue("font_scale",font_scale,
			"Set the scaling of the font used for labeling "
			"channels. A value of 1.0 means the labels should "
			"fill the channel diameter. Smaller values mean the "
			"labels will be smaller in the channel.");

  options.findWithValue("colormap",colormap,
			"Set the color map (palette) to use when "
			"\"-valuemode=color\". Valid options are "
			"\"jet\", \"hot\" and \"blue_red_yellow\".");

  options.findWithValue("colors",colors,
			"Directly set the colors to use. Colors should be "
			"given as a comma separated list of strings, each "
			"of which should be a color as understood by QT. This "
			"option overrides the \"colormap\" option.");

  options.findWithValue("png_file", png_file,
			"Name of output PNG file, empty string suppresses PNG "
			"output.");

  options.findWithValue("png_size", png_size,
			"Set width and height of PNG image in pixels.");

  options.findWithValue("ps_file", ps_file,
			"Name of output PS file, empty string suppresses PS "
			"output.");

#if ( QT_VERSION >= 0x040100 )
  options.findWithValue("pdf_file", pdf_file,
			"Name of output PDF file, empty string suppresses PDF "
			"output.");
#endif

  options.findWithValue("value_column", value_column,
			"Select the column to display (counting from zero). "
			"This option is overridden by \"-columns\".");

  options.findWithValue("value_default", value_default,
			"Select the default value in no data is given");

  columns.insert(columns.end(), value_column, "");
  columns.push_back("value");

  options.findWithValue("columns", columns,
			"Define contents of columns in input data. The "
			"columns should be specified as a comma separated "
			"list containing some or all of the following items "
			"in any order: \"scope\", \"chan\", \"value\", "
			"\"type\", \"suppress\", \"highlight\", \"no_draw\" "
			"and \"no_data\".");

  if(options.find("log_scale", "Set log scale.") != VSOptions::FS_NOT_FOUND)
    log_scale = true;

  if(options.find("subtract_mean", "Subtract mean from each tube value.")
     !=VSOptions::FS_NOT_FOUND)
    subtract_mean = true;

  options.findWithValue("limits", limits, 
			"Fix limits on scale. Limits should be given in the "
			"form: \"min/max\". Automatic scaling is done if "
			"min=max.");

  if(options.find("clip_value", "Clip data value to limits.")
     !=VSOptions::FS_NOT_FOUND)
    clip_data = true;

  std::string progname(*my_argv);
  my_argv++,my_argc--;

  if(!options.assertNoOptions())
    {
      std::cerr << progname << ": unknown options: ";
      for(int i=1;i<my_argc;i++)
        if(*(my_argv[i])=='-') std::cerr << ' ' << my_argv[i];
      std::cerr << std::endl;
      usage(progname, options, std::cerr);
      exit(EXIT_FAILURE);
    }

  if(print_usage)
    {
      usage(progname, options, std::cout);
      exit(EXIT_SUCCESS);
    }

  if(my_argc>1)
    {
      std::cerr << progname
		<< ": too many argumenets given (" << my_argc << ')'
		<< std::endl << std::endl;
      usage(progname, options, std::cerr);
      exit(EXIT_FAILURE);
    }

  std::string filename;
  if(my_argc)
    {
      filename = *my_argv;
      my_argv++,my_argc--;
    }

  VSimpleCameraLayout* camera;
  
  camera_type = toUpper(camera_type);
  if(camera_type == "STANDARD")
    {
      VSimpleCameraFactory::WhippleCamera wcamera = 
	VSimpleCameraFactory::VC499BIGTUBES;
      whipple_camera = toUpper(whipple_camera);
      if(whipple_camera == "WC109")wcamera=VSimpleCameraFactory::WC109;
      else if(whipple_camera == "WC151")wcamera=VSimpleCameraFactory::WC151;
      else if(whipple_camera == "WC331")wcamera=VSimpleCameraFactory::WC331;
      else if(whipple_camera == "WC490")wcamera=VSimpleCameraFactory::WC490;
      else if(whipple_camera == "VC499")wcamera=VSimpleCameraFactory::VC499;
      else if(whipple_camera == "VC499BIGTUBES")
	wcamera=VSimpleCameraFactory::VC499BIGTUBES;
      else if(whipple_camera == "VC500")wcamera=VSimpleCameraFactory::VC500;
      else if(whipple_camera == "VC500BIGTUBES")
	wcamera=VSimpleCameraFactory::VC500BIGTUBES;
      else
	{
	  std::cerr << progname << ": unknown standard camera: " 
		    << whipple_camera << std::endl;
	  std::cerr << std::endl;
	  usage(progname, options, std::cerr);
	  exit(EXIT_FAILURE);      
	}
      camera = VSimpleCameraFactory::getInstance()->
	getWhippleCamera(wcamera);
    }
  else if(camera_type == "ROUNDED")
    {
      camera = VSimpleCameraFactory::getInstance()->
	getRoundedHexCamera(camera_diameter/2.0, tube_diameter/2.0);
    }
  else if(camera_type == "HEX")
    {
      camera = VSimpleCameraFactory::getInstance()->
	getHexCamera(hex_rings, tube_diameter/2.0);
    }
  else
    {
      std::cerr << progname << ": unknown camera type: " << camera_type
		<< std::endl;
      std::cerr << std::endl;
      usage(progname, options, std::cerr);
      exit(EXIT_FAILURE);      
    }

  VQMultiCamera::camera_list cameras;
  cameras.resize(ncam,camera);

  VQMCArrayData array_data(ncam);
  std::vector<std::vector<double> > raw_values(ncam);
  for(unsigned icam=0;icam!=ncam;icam++)
    {
      unsigned nchan = cameras[icam]->numChannels();
      array_data.fCameraData[icam].fChannelData.resize(nchan, 
						       VQMCChannelData(true));
      raw_values[icam].resize(nchan);
    }

  std::istream* datastream;
  std::ifstream filestream;

  if(filename.empty())
    {
      datastream = &std::cin;
      filename = "<stdin>";
    }
  else
    {
      filestream.open(filename.c_str());
      datastream = &filestream;
    }

  unsigned icam = 0;
  unsigned ichan = 0;
  std::string line;
  VSLineTokenizer tokenizer;
  VSTokenList tokens;
  unsigned iline=0;
  while(getline(*datastream,line))
    {
      iline++;
      std::string line_copy = line;
      tokenizer.tokenize(line_copy, tokens);
      if(tokens.size() == 0)continue;
      if(tokens.size() < columns.size())
	{
	  std::cerr << progname << ": " << filename << ": line " << iline 
		    << ": " << tokens.size() << "columns found (" 
		    << columns.size() << " expected)" << std::endl
		    << "Line: " << line << std::endl;
	  continue;
	}

      double value = value_default;
      unsigned type = 0;
      bool suppress = false;
      bool highlight = false;
      bool no_draw = false;
      bool no_data = false;

      bool found_chan = false;

      for(unsigned icolumn = 0; icolumn<columns.size(); icolumn++)
	{
	  if(columns[icolumn]=="")
	    continue;
	  else if(columns[icolumn]=="scope")
	    tokens[icolumn].to(icam);
	  else if(columns[icolumn]=="chan")
	    tokens[icolumn].to(ichan),found_chan=true;
	  else if(columns[icolumn]=="value")
	    tokens[icolumn].to(value);
	  else if(columns[icolumn]=="type")
	    tokens[icolumn].to(type);
	  else if(columns[icolumn]=="suppress")
	    tokens[icolumn].to(suppress);
	  else if(columns[icolumn]=="highlight")
	    tokens[icolumn].to(highlight);
	  else if(columns[icolumn]=="no_draw")
	    tokens[icolumn].to(no_draw);
	  else if(columns[icolumn]=="no_data")
	    tokens[icolumn].to(no_data);
	  else
	    {
	      std::cerr << progname << ": unknown column: " << columns[icolumn]
			<< std::endl;
	      std::cerr << std::endl;
	      usage(progname, options, std::cerr);
	      exit(EXIT_FAILURE);      
	    }
	}

      if((icam<ncam)&&(ichan<array_data.fCameraData[icam].fChannelData.size()))
	{
	  raw_values[icam][ichan]=value;
	  if(log_scale)
	    {
	      if(value>0)value=log10(value);
	      else value=0,no_data=true;
	    }

	  VQMCChannelData& 
	    data(array_data.fCameraData[icam].fChannelData[ichan]);
	  data.fValue = value;
	  data.fType = type;
	  data.fSuppress = suppress;
	  data.fHighlight = highlight;
	  data.fNoDraw = no_draw;
	  data.fNoData = no_data;
	}
      else
	{
	  std::cerr << progname << ": " << filename << ": line " << iline 
		    << ": unknown scope/channel: " << icam << '/' << ichan 
		    << std::endl
		    << "Line: " << line << std::endl;
	}

      if(found_chan)continue;

      ichan++;
      if(ichan >= cameras[icam]->numChannels())
	{
	  ichan=0;
	  icam++;
	  if(icam>=ncam)break;
	}
    }
  
  if(filename != "<stdin>")filestream.close();

  bool found_data = false;
  unsigned min_type = 0;
  unsigned max_type = 0;
  for(unsigned icam=0; icam<ncam; icam++)
    {
      unsigned nchan = array_data.fCameraData[icam].fChannelData.size();
      for(unsigned ichan=0; ichan<nchan; ichan++)
	{
	  VQMCChannelData& 
	      data(array_data.fCameraData[icam].fChannelData[ichan]);
	  if(!data.fNoData)
	    {
	      if((!found_data)||(data.fType<min_type))min_type=data.fType;
	      if((!found_data)||(data.fType>max_type))max_type=data.fType;
	      found_data=true;
	    }
	}
    }  

  if(subtract_mean)
    for(unsigned icam=0; icam<ncam; icam++)
      {
	double sum = 0;
	unsigned count = 0;
	unsigned nchan = array_data.fCameraData[icam].fChannelData.size();
	for(unsigned ichan=0; ichan<nchan; ichan++)
	  {
	    VQMCChannelData& 
	      data(array_data.fCameraData[icam].fChannelData[ichan]);
	    if(!data.fNoData)sum+=data.fValue, count++;
	  }
	double mean = sum/double(count);
	for(unsigned ichan=0; ichan<nchan; ichan++)
	  {
	    VQMCChannelData& 
	      data(array_data.fCameraData[icam].fChannelData[ichan]);
	    if(!data.fNoData)data.fValue-=mean;
	  }
      }

  VQMCSimpleChannelRendererBase::ValueMode VM;
  value_mode = toLower(value_mode);
  if(value_mode == "area")
    VM=VQMCSimpleChannelRendererBase::VM_AREA;
  else if(value_mode == "radius")
    VM=VQMCSimpleChannelRendererBase::VM_RADIUS;
  else if(value_mode == "color")
    VM=VQMCSimpleChannelRendererBase::VM_COLOR;
  else 
    {
      std::cerr << progname << ": unknown value mode: " << value_mode
		<< std::endl;
      std::cerr << std::endl;
      usage(progname, options, std::cerr);
      exit(EXIT_FAILURE);      
    }

  VQMCSimpleChannelRendererBase::ColorMode CM;
  color_mode = toLower(color_mode);
  if(color_mode == "scope")
    CM=VQMCSimpleChannelRendererBase::CM_BY_CAMERA;
  else if(color_mode == "type")
    CM=VQMCSimpleChannelRendererBase::CM_BY_TYPE;
  else 
    {
      std::cerr << progname << ": unknown color mode: " << color_mode
		<< std::endl;
      std::cerr << std::endl;
      usage(progname, options, std::cerr);
      exit(EXIT_FAILURE);      
    }

  std::vector<QColor> the_colormap;
  colormap = toLower(colormap);
  the_colormap = VQColormapFactory::getColormap(colormap);
  if(the_colormap.empty())
    {
      std::cerr << progname << ": unknown colormap: " << colormap
		<< std::endl;
      std::cerr << std::endl;
      usage(progname, options, std::cerr);
      exit(EXIT_FAILURE);      
    }

  std::vector<QColor> COL;
  if(!colors.empty())
    {
      for(std::vector<std::string>::const_iterator icol = colors.begin();
	  icol != colors.end(); icol++)
	{
	  QColor col(icol->c_str());
	  if(col.isValid())COL.push_back(col);
	  else std::cerr << progname << ": warning: invalid color: "
			 << *icol << std::endl;
	}
    }
  else if(VM == VQMCSimpleChannelRendererBase::VM_COLOR)
    {
      COL = the_colormap;
    }
  else
    {
      unsigned ncol=ncam;
      if(CM==VQMCSimpleChannelRendererBase::CM_BY_TYPE)ncol=max_type+1;
      
      for(unsigned icol=0; icol<ncol; icol++)
	{
	  double y = double(icol)/double(ncol);
	  if(y>1)y=1;
	  else if(y<0)y=0;
	  double z = y*double(the_colormap.size()-1);
	  unsigned index = unsigned(floor(z));
	  double d = fmod(z,1);
	  QColor color = the_colormap[index];
	  if(index != the_colormap.size()-1)
	    blendColors(1-d,color,the_colormap[index+1]);
	  COL.push_back(color);
	}
    }

  bool DL = !no_draw_outline;

  VQMCChannelRenderer* channel_renderer = 0;
  draw_type = toUpper(draw_type);
  
  if(draw_type == "CIRCLE")
    channel_renderer =
      new VQMCSimpleChannelRenderer<VQMCCircleDraw>(COL, CM, VM, DL);
  else if(draw_type == "SQUARE")
    channel_renderer =
      new VQMCSimpleChannelRenderer<VQMCSquareDraw>(COL, CM, VM, DL);
  else if(draw_type == "SMALLSQUARE")
    channel_renderer =
      new VQMCSimpleChannelRenderer<VQMCSmallSquareDraw>(COL, CM, VM, DL);
  else if(draw_type == "HEX")
    channel_renderer =
      new VQMCSimpleChannelRenderer<VQMCHexDraw>(COL, CM, VM, DL);
  else if(draw_type == "SMALLHEX")
    channel_renderer =
      new VQMCSimpleChannelRenderer<VQMCSmallHexDraw>(COL, CM, VM, DL);
  else if(draw_type == "STAR")
    channel_renderer =
      new VQMCSimpleChannelRenderer<VQMCStarDraw>(COL, CM, VM, DL);
  else if(draw_type == "PIE")
    channel_renderer =
      new VQMCSimpleChannelRenderer<VQMCPieDraw>(COL, CM, VM, DL);
  else
    {
      std::cerr << progname << ": unknown tube draw type: " << draw_type
		<< std::endl;
      std::cerr << std::endl;
      usage(progname, options, std::cerr);
      exit(EXIT_FAILURE);      
    }

  VQMultiCamera* widget = 0;
  VQMultiCameraRenderer* renderer;

  if(!batch)
    {
      widget = new VQMultiCamera(channel_renderer,cameras);
      renderer = widget;
      widget->setFastRender(draw_fast);
    }
  else
    {
      renderer = new VQMultiCameraRenderer(channel_renderer,cameras);
    }

  //renderer->setChannelRenderMode(VQMultiCamera::CRM_MAX_VALUE_PER_CHANNEL);
  if(limits.first != limits.second)
    renderer->setAllCameraLimits(VQMCCameraDataLimits(limits.first,
						      limits.second));
  renderer->setClipToLimits(clip_data);
  renderer->setUseCommonLimits(true);

  renderer->setLabelChannels(label_channels);
  renderer->setFontFamily(font_family);
  renderer->setFontScale(font_scale);
  
  unsigned nx_cam = unsigned(ceil(sqrt(double(ncam))));
  unsigned ny_cam = unsigned(ceil(double(ncam)/double(nx_cam)));

  if(!camera_overlap)
    {
      double max_x_extent = 0;
      double max_y_extent = 0;
      for(unsigned icam=0;icam<ncam;icam++)
	{
	  double ex = cameras[icam]->extentRight()-cameras[icam]->extentLeft();
	  double ey = cameras[icam]->extentTop()-cameras[icam]->extentBottom();
	  if(max_x_extent < ex)max_x_extent=ex;
	  if(max_y_extent < ey)max_y_extent=ey;
	}
      max_x_extent *= 1.02;
      max_y_extent *= 1.02;
      
      for(unsigned icam=0;icam<ncam;icam++)
	{
	  unsigned ix = icam%nx_cam;
	  unsigned iy = icam/nx_cam;
	  double x = 
	    (cameras[icam]->extentRight()+cameras[icam]->extentLeft())/2;
	  double y = 
	    (cameras[icam]->extentTop()+cameras[icam]->extentBottom())/2;
	  x += (double(ix)-double(nx_cam-1)/2)*max_x_extent;
	  y += (-double(iy)+double(ny_cam-1)/2)*max_y_extent;
	  renderer->setCameraOffsetXDeg(icam, x);
	  renderer->setCameraOffsetYDeg(icam, y);
	}
    }

  if(!batch)
    {
      widget->setArrayData(array_data);
      widget->show();
      app.connect(&app,SIGNAL(lastWindowClosed()),&app,SLOT(quit()));
      if(!no_pause)
	app.exec();
    }

  if(!ps_file.empty())
    {
      QPrinter* output = new QPrinter(QPrinter::HighResolution);
      output->setOutputFileName(ps_file.c_str());
      QPainter painter(output);
      renderer->configurePainter(painter, false);
      renderer->renderCamera(painter,array_data, false);
      painter.end();
      delete output;
    }

#if ( QT_VERSION >= 0x040100 )
  if(!pdf_file.empty())
    {
      QPrinter* output = new QPrinter(QPrinter::HighResolution);
      output->setOutputFormat(QPrinter::PdfFormat);
      output->setOutputFileName(pdf_file.c_str());
      QPainter painter(output);
      renderer->configurePainter(painter, false);
      renderer->renderCamera(painter,array_data, false);
      painter.end();
      delete output;
    }
#endif

  if(!png_file.empty())
    {
      QImage* output = 
	new QImage(png_size,png_size,QImage::Format_ARGB32_Premultiplied);
      QPainter painter(output);
      painter.setBackground(QBrush(Qt::white));
      painter.fillRect(output->rect(),painter.background());
      renderer->configurePainter(painter, false);
      renderer->renderCamera(painter,array_data, false);
      painter.end();
      output->save(png_file.c_str(),"PNG");
      delete output;
    }

  //delete[] my_argv;
  delete renderer;
  delete channel_renderer;
}
