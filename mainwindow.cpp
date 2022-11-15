#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "filter.h"
#include "qcustomplot.h"


typedef struct
{
  unsigned numStages;         /**< number of 2nd order stages in the filter.  Overall order is 2*numStages. */
  float *pState;          /**< Points to the array of state coefficients.  The array is of length 4*numStages. */
  float *pCoeffs;         /**< Points to the array of coefficients.  The array is of length 5*numStages. */
} arm_biquad_casd_df1_inst_f32;

#define N 4

//float coaf[5*N] = {
//    -2.192418628271423750, 0.000000000000000000, 2.192418628271423750, +1.987056911885394460, -0.987106041330942929,
//    -0.101505574010362101, 0.000000000000000000, 0.101505574010362101, +0.154870658096936065, -0.045701475093464263,
//    -0.823371366424797357, 0.000000000000000000, 0.823371366424797357, +1.994650922359786360, -0.994699676722548332,
//    -0.375225408105429370, 0.000000000000000000, 0.375225408105429370, +0.207072761974852609, -0.453302856232669915,
//};

//float coaf[5*N]={
//     1.12796421E+00, 2.17021499E+00, 1.12796421E+00, 1.98860953E+00, -9.88675700E-01,
//     1.12796421E+00, -2.25592822E+00, 1.12796421E+00, 1.06500021E+00, -3.83519397E-01,
//     1.26058951E-01, 2.01750553E-01, 1.26058951E-01, 7.80095015E-01, -7.43266042E-01,
//     1.26058951E-01, -2.52117780E-01, 1.26058951E-01, 1.99841009E+00, -9.98432752E-01,};

float coaf[5*N]={ 1.03417331E+00, 1.96991496E+00, 1.03417331E+00, 1.98862045E+00, -9.88686679E-01,
                  1.03417331E+00, -2.06834645E+00, 1.03417331E+00, 1.15577315E+00, -4.22240356E-01,
                  1.17588826E-01, 1.77532018E-01, 1.17588826E-01, 9.38136820E-01, -7.55808691E-01,
                  1.17588826E-01, -2.35177540E-01, 1.17588826E-01, 1.99841143E+00, -9.98434093E-01,};

float state[N*4];
arm_biquad_casd_df1_inst_f32 S;


void arm_biquad_cascade_df1_init_f32
(arm_biquad_casd_df1_inst_f32 * S, unsigned char numStages, float * pCoeffs, float * pState)
{
  /* Assign filter stages */
  S->numStages = numStages;

  /* Assign coefficient pointer */
  S->pCoeffs = pCoeffs;

  /* Clear state buffer and size is always 4 * numStages */
  memset(pState, 0, (4u * (unsigned) numStages) * sizeof(float));

  /* Assign state pointer */
  S->pState = pState;
}

void arm_biquad_cascade_df1_f32
(const arm_biquad_casd_df1_inst_f32 * S, float * pSrc, float * pDst, unsigned blockSize)
{
  float *pIn = pSrc;                         /*  source pointer            */
  float *pOut = pDst;                        /*  destination pointer       */
  float *pState = S->pState;                 /*  pState pointer            */
  float *pCoeffs = S->pCoeffs;               /*  coefficient pointer       */
  float acc;                                 /*  Simulates the accumulator */
  float b0, b1, b2, a1, a2;                  /*  Filter coefficients       */
  float Xn1, Xn2, Yn1, Yn2;                  /*  Filter pState variables   */
  float Xn;                                  /*  temporary input           */
  unsigned sample, stage = S->numStages;         /*  loop counters             */

  /* Run the below code for Cortex-M4 and Cortex-M3 */
  do
  {
    /* Reading the coefficients */
    b0 = *pCoeffs++;
    b1 = *pCoeffs++;
    b2 = *pCoeffs++;
    a1 = *pCoeffs++;
    a2 = *pCoeffs++;

    /* Reading the pState values */
    Xn1 = pState[0];
    Xn2 = pState[1];
    Yn1 = pState[2];
    Yn2 = pState[3];

    /* Apply loop unrolling and compute 4 output values simultaneously. */
    /*      The variable acc hold output values that are being computed:
     *
     *    acc =  b0 * x[n] + b1 * x[n-1] + b2 * x[n-2] + a1 * y[n-1]   + a2 * y[n-2]
     *    acc =  b0 * x[n] + b1 * x[n-1] + b2 * x[n-2] + a1 * y[n-1]   + a2 * y[n-2]
     *    acc =  b0 * x[n] + b1 * x[n-1] + b2 * x[n-2] + a1 * y[n-1]   + a2 * y[n-2]
     *    acc =  b0 * x[n] + b1 * x[n-1] + b2 * x[n-2] + a1 * y[n-1]   + a2 * y[n-2]
     */

    sample = blockSize >> 2u;

    /* First part of the processing with loop unrolling.  Compute 4 outputs at a time.
     ** a second loop below computes the remaining 1 to 3 samples. */
    while(sample > 0u)
    {
      /* Read the first input */
      Xn = *pIn++;

      /* acc =  b0 * x[n] + b1 * x[n-1] + b2 * x[n-2] + a1 * y[n-1] + a2 * y[n-2] */
      Yn2 = (b0 * Xn) + (b1 * Xn1) + (b2 * Xn2) + (a1 * Yn1) + (a2 * Yn2);

      /* Store the result in the accumulator in the destination buffer. */
      *pOut++ = Yn2;

      /* Every time after the output is computed state should be updated. */
      /* The states should be updated as:  */
      /* Xn2 = Xn1    */
      /* Xn1 = Xn     */
      /* Yn2 = Yn1    */
      /* Yn1 = acc   */

      /* Read the second input */
      Xn2 = *pIn++;

      /* acc =  b0 * x[n] + b1 * x[n-1] + b2 * x[n-2] + a1 * y[n-1] + a2 * y[n-2] */
      Yn1 = (b0 * Xn2) + (b1 * Xn) + (b2 * Xn1) + (a1 * Yn2) + (a2 * Yn1);

      /* Store the result in the accumulator in the destination buffer. */
      *pOut++ = Yn1;

      /* Every time after the output is computed state should be updated. */
      /* The states should be updated as:  */
      /* Xn2 = Xn1    */
      /* Xn1 = Xn     */
      /* Yn2 = Yn1    */
      /* Yn1 = acc   */

      /* Read the third input */
      Xn1 = *pIn++;

      /* acc =  b0 * x[n] + b1 * x[n-1] + b2 * x[n-2] + a1 * y[n-1] + a2 * y[n-2] */
      Yn2 = (b0 * Xn1) + (b1 * Xn2) + (b2 * Xn) + (a1 * Yn1) + (a2 * Yn2);

      /* Store the result in the accumulator in the destination buffer. */
      *pOut++ = Yn2;

      /* Every time after the output is computed state should be updated. */
      /* The states should be updated as: */
      /* Xn2 = Xn1    */
      /* Xn1 = Xn     */
      /* Yn2 = Yn1    */
      /* Yn1 = acc   */

      /* Read the forth input */
      Xn = *pIn++;

      /* acc =  b0 * x[n] + b1 * x[n-1] + b2 * x[n-2] + a1 * y[n-1] + a2 * y[n-2] */
      Yn1 = (b0 * Xn) + (b1 * Xn1) + (b2 * Xn2) + (a1 * Yn2) + (a2 * Yn1);

      /* Store the result in the accumulator in the destination buffer. */
      *pOut++ = Yn1;

      /* Every time after the output is computed state should be updated. */
      /* The states should be updated as:  */
      /* Xn2 = Xn1    */
      /* Xn1 = Xn     */
      /* Yn2 = Yn1    */
      /* Yn1 = acc   */
      Xn2 = Xn1;
      Xn1 = Xn;

      /* decrement the loop counter */
      sample--;

    }

    /* If the blockSize is not a multiple of 4, compute any remaining output samples here.
     ** No loop unrolling is used. */
    sample = blockSize & 0x3u;

    while(sample > 0u)
    {
      /* Read the input */
      Xn = *pIn++;

      /* acc =  b0 * x[n] + b1 * x[n-1] + b2 * x[n-2] + a1 * y[n-1] + a2 * y[n-2] */
      acc = (b0 * Xn) + (b1 * Xn1) + (b2 * Xn2) + (a1 * Yn1) + (a2 * Yn2);

      /* Store the result in the accumulator in the destination buffer. */
      *pOut++ = acc;

      /* Every time after the output is computed state should be updated. */
      /* The states should be updated as:    */
      /* Xn2 = Xn1    */
      /* Xn1 = Xn     */
      /* Yn2 = Yn1    */
      /* Yn1 = acc   */
      Xn2 = Xn1;
      Xn1 = Xn;
      Yn2 = Yn1;
      Yn1 = acc;

      /* decrement the loop counter */
      sample--;

    }

    /*  Store the updated state variables back into the pState array */
    *pState++ = Xn1;
    *pState++ = Xn2;
    *pState++ = Yn1;
    *pState++ = Yn2;

    /*  The first stage goes from the input buffer to the output buffer. */
    /*  Subsequent numStages  occur in-place in the output buffer */
    pIn = pDst;

    /* Reset the output pointer */
    pOut = pDst;

    /* decrement the loop counter */
    stage--;

  } while(stage > 0u);
}

template<typename T>
T getMin(QVector<T> data, int dataSize)
{
    int smallest_element_raw = data[0]; //let, first element is the smallest one
    for(int i = 1; i < dataSize; i++)  //start iterating from the second element
    {
        if(data[i] < smallest_element_raw)
        {
           smallest_element_raw = data[i];
        }
    }
    return smallest_element_raw;
}

template<typename T>
T getMax(QVector<T> data, int dataSize)
{
    int largest_element_raw = data[0]; //also let, first element is the biggest one
    for(int i = 1; i < dataSize; i++)  //start iterating from the second element
    {
        if(data[i] > largest_element_raw)
        {
           largest_element_raw = data[i];
        }
    }
    return largest_element_raw;
}


void arm_biquad_cascade_df1_init_f32
(arm_biquad_casd_df1_inst_f32 * S, unsigned char numStages, float * pCoeffs, float * pState);

void arm_biquad_cascade_df1_f32
(const arm_biquad_casd_df1_inst_f32 * S, float * pSrc, float * pDst, unsigned blockSize);

void sweep(double f_start, double f_end, double interval, int n_steps) {
    for (int i = 0; i < n_steps; ++i) {
        double delta = i / (float)n_steps;
        double t = interval * delta;
        double phase = 2 * M_PI * t * (f_start + (f_end - f_start) * delta / 2);
        while (phase > 2 * M_PI) phase -= 2 * M_PI; // optional
        printf("%f %f %f", t, phase * 180 / M_PI, 3 * qSin(phase));
    }
}

const double TwoPi = 6.283185307179586;
void FFTAnalysis(double *AVal, double *FTvl, int Nvl, int Nft) {
  int i, j, n, m, Mmax, Istp;
  double Tmpr, Tmpi, Wtmp, Theta;
  double Wpr, Wpi, Wr, Wi;
  double *Tmvl;

  n = Nvl * 2; Tmvl = new double[n];

  for (i = 0; i < n; i+=2) {
   Tmvl[i] = 0;
   Tmvl[i+1] = AVal[i/2];
  }

  i = 1; j = 1;
  while (i < n) {
    if (j > i) {
      Tmpr = Tmvl[i]; Tmvl[i] = Tmvl[j]; Tmvl[j] = Tmpr;
      Tmpr = Tmvl[i+1]; Tmvl[i+1] = Tmvl[j+1]; Tmvl[j+1] = Tmpr;
    }
    i = i + 2; m = Nvl;
    while ((m >= 2) && (j > m)) {
      j = j - m; m = m >> 1;
    }
    j = j + m;
  }

  Mmax = 2;
  while (n > Mmax) {
    Theta = -TwoPi / Mmax; Wpi = sin(Theta);
    Wtmp = sin(Theta / 2); Wpr = Wtmp * Wtmp * 2;
    Istp = Mmax * 2; Wr = 1; Wi = 0; m = 1;

    while (m < Mmax) {
      i = m; m = m + 2; Tmpr = Wr; Tmpi = Wi;
      Wr = Wr - Tmpr * Wpr - Tmpi * Wpi;
      Wi = Wi + Tmpr * Wpi - Tmpi * Wpr;

      while (i < n) {
        j = i + Mmax;
        Tmpr = Wr * Tmvl[j] - Wi * Tmvl[j-1];
        Tmpi = Wi * Tmvl[j] + Wr * Tmvl[j-1];

        Tmvl[j] = Tmvl[i] - Tmpr; Tmvl[j-1] = Tmvl[i-1] - Tmpi;
        Tmvl[i] = Tmvl[i] + Tmpr; Tmvl[i-1] = Tmvl[i-1] + Tmpi;
        i = i + Istp;
      }
    }

    Mmax = Istp;
  }

  for (i = 0; i < Nft; i++) {
    j = i * 2; FTvl[i] = 2*sqrt(pow(Tmvl[j],2) + pow(Tmvl[j+1],2))/Nvl;
  }

  delete []Tmvl;
}

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    arm_biquad_cascade_df1_init_f32(&S, N,(float*)coaf, state);
    on_pushButtonFilter_clicked();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::plotGrapics()
{
    int SAMPLE_COUNT = ui->lineEdit_sampleRate->text().toInt();
    ui->horizontalSlider->setMaximum(SAMPLE_COUNT);

    float* bufferDestination = new float[SAMPLE_COUNT];
    float* bufferSource = new float[SAMPLE_COUNT];

    const int sineFreq = ui->lineEdit_sineFreq->text().toInt();

    unsigned spectre_size = log(SAMPLE_COUNT)/log(2);
        spectre_size = pow(2, ++spectre_size);

    QVector<double> x_all, y_raw, y_filtred, y_magnitude;

    double* Input_spectre = new double[spectre_size];
    double* Output_spectre = new double[spectre_size];
    for(unsigned i = 0; i < spectre_size; i++)
    {
        Input_spectre[0] = 0;
        Output_spectre[0] = 0;
    }

    for(int i = 0; i < SAMPLE_COUNT; i++)
    {
        bufferSource[i] = qSin(2 * M_PI * i * sineFreq / SAMPLE_COUNT);
    }

    arm_biquad_cascade_df1_f32(&S, bufferSource, bufferDestination, SAMPLE_COUNT);

    for(int i = 0; i < SAMPLE_COUNT; i++)
    {
        x_all.push_back(i);

        y_raw.push_back(bufferSource[i]);
        y_filtred.push_back(bufferDestination[i]);

        Input_spectre[i] = bufferDestination[i];
    }

    FFTAnalysis(Input_spectre, Output_spectre, spectre_size, spectre_size);

    for(int i = 0; i < SAMPLE_COUNT / 2; i++)
    {
        y_magnitude.push_back(Output_spectre[i]);
    }

    delete[] Input_spectre;
    delete[] Output_spectre;
    delete[] bufferSource;
    delete[] bufferDestination;

    QFont legendFont = font();
    legendFont.setPointSize(8);

    /* --------------------RAW SINE-------------------- */
    ui->canvas_raw->clearGraphs();
    ui->canvas_raw->xAxis->setRange(0, SAMPLE_COUNT/ sineFreq * 2);
    ui->canvas_raw->yAxis->setRange(-1.5, 1.5);

    ui->canvas_raw->xAxis->setAutoTickCount(15);
    ui->canvas_raw->yAxis->setAutoTickCount(10);

    ui->canvas_raw->legend->clear();
    ui->canvas_raw->legend->setVisible(true);
    ui->canvas_raw->legend->setFont(legendFont);
    ui->canvas_raw->legend->setBrush(QBrush(QColor(255,255,255,230)));

    ui->canvas_raw->xAxis->setLabel("Units");
    ui->canvas_raw->yAxis->setLabel("Units");

    ui->canvas_raw->addGraph();
    ui->canvas_raw->graph(0)->setData(x_all, y_raw);
    ui->canvas_raw->graph(0)->setName(QString("Sine before filtration"));
    ui->canvas_raw->replot();

    ui->canvas_raw->setInteraction(QCP::iRangeDrag, true);
    ui->canvas_raw->setInteraction(QCP::iRangeZoom, true);

    /* --------------------FILTRED SINE-------------------- */
    ui->canvas_filtred->clearGraphs();
    ui->canvas_filtred->xAxis->setRange(0, SAMPLE_COUNT / sineFreq * 2);
    ui->canvas_filtred->yAxis->setRange(-1.5, 1.5);

    ui->canvas_filtred->xAxis->setAutoTickCount(15);
    ui->canvas_filtred->yAxis->setAutoTickCount(10);

    ui->canvas_filtred->legend->clear();
    ui->canvas_filtred->legend->setVisible(true);
    ui->canvas_filtred->legend->setFont(legendFont);
    ui->canvas_filtred->legend->setBrush(QBrush(QColor(255,255,255,230)));

    ui->canvas_filtred->xAxis->setLabel("Units");
    ui->canvas_filtred->yAxis->setLabel("Units");

    ui->canvas_filtred->addGraph();
    ui->canvas_filtred->graph(0)->setData(x_all, y_filtred);
    ui->canvas_filtred->graph(0)->setName(QString("Sine after filtration"));
    ui->canvas_filtred->replot();

    ui->canvas_filtred->setInteraction(QCP::iRangeDrag, true);
    ui->canvas_filtred->setInteraction(QCP::iRangeZoom, true);

    /* --------------------IMPULSE-------------------- */
    ui->canvas_impulse->clearGraphs();
    ui->canvas_impulse->xAxis->setRange(-100, SAMPLE_COUNT / 2 + 100);
    ui->canvas_impulse->yAxis->setRange(-0.5, 1);

    ui->canvas_impulse->xAxis->setAutoTickCount(15);
    ui->canvas_impulse->yAxis->setAutoTickCount(20);

    ui->canvas_impulse->legend->clear();
    ui->canvas_impulse->legend->setVisible(true);
    ui->canvas_impulse->legend->setFont(legendFont);
    ui->canvas_impulse->legend->setBrush(QBrush(QColor(255,255,255,230)));

    ui->canvas_impulse->xAxis->setLabel("Hz");
    ui->canvas_impulse->yAxis->setLabel("Units");

    ui->canvas_impulse->addGraph();
    ui->canvas_impulse->graph(0)->setData(x_all, y_magnitude);
    ui->canvas_impulse->graph(0)->setName(QString("Spectrum of sine"));
    ui->canvas_impulse->replot();

    ui->canvas_impulse->setInteraction(QCP::iRangeDrag, true);
    ui->canvas_impulse->setInteraction(QCP::iRangeZoom, true);
}

void MainWindow::on_pushButtonFilter_clicked()
{
    plotGrapics();
}

void MainWindow::on_horizontalSlider_valueChanged(int value)
{
   ui->lineEdit_sineFreq->setText(QString("%1").arg(QString::number(value)));
   plotGrapics();
}
