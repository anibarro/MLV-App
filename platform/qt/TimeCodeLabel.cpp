/*!
 * \file TimeCodeLabel.cpp
 * \author masc4ii
 * \copyright 2018
 * \brief Paint TimeCode Image
 */

#include "TimeCodeLabel.h"
#include <QRect>
#include <QPainter>
#include <QLinearGradient>
#include <QFontDatabase>
#include <QString>

//Constructor
TimeCodeLabel::TimeCodeLabel()
{
    m_tcImage = new QImage( 400, 60, QImage::Format_RGB888 );

    QFontDatabase::addApplicationFont( ":/Fonts/Fonts/DSEG7Modern-Regular.ttf" );
}

//Destructor
TimeCodeLabel::~TimeCodeLabel()
{
    delete m_tcImage;
}

//Get the image for frameNumber at clipFps
QImage TimeCodeLabel::getTimeCodeLabel(uint32_t frameNumber, float clipFps)
{
    //Background
    m_tcImage->fill( QColor( 20, 20, 20, 255 ) );

    //Boarder
    QRect rect( 0, 0, 2, 60 );
    QPainter painterTc( m_tcImage );
    QLinearGradient gradient( rect.topLeft(), (rect.topLeft() + rect.bottomLeft() ) / 2 ); // diagonal gradient from top-left to bottom-right
    gradient.setColorAt( 0, QColor( 30, 30, 30, 255 ) );
    gradient.setColorAt( 1, QColor( 90, 90, 90, 255 ) );
    painterTc.fillRect(rect, gradient);

    QRect rect2( 398, 0, 2, 60 );
    painterTc.fillRect(rect2, gradient);

    QRect rect3( 0, 0, 400, 2 );
    QLinearGradient gradient3( rect3.topLeft(), rect3.bottomLeft() );
    gradient3.setColorAt( 0, QColor( 90, 90, 90, 255 ) );
    gradient3.setColorAt( 1, QColor( 30, 30, 30, 255 ) );
    painterTc.fillRect(rect3, gradient3);

    QRect rect4( 0, 58, 400, 2 );
    QLinearGradient gradient4( rect4.topLeft(), rect4.bottomLeft() );
    gradient4.setColorAt( 0, QColor( 90, 90, 90, 255 ) );
    gradient4.setColorAt( 1, QColor( 60, 60, 60, 255 ) );
    painterTc.fillRect(rect4, gradient4);

    //Text
    int frame = (int) (frameNumber % (int)clipFps);
    frameNumber /= (int)clipFps;
    int seconds = (int) (frameNumber % 60);
    frameNumber /= 60;
    int minutes = (int) (frameNumber % 60);
    frameNumber /= 60;
    int hours = (int) (frameNumber);

    //Font selection
#ifdef Q_OS_MAC
    QFont font = QFont("DSEG7 Modern", 32, 1);
#else
    QFont font = QFont("DSEG7 Modern", 24, 1);
#endif

    //Shadow
    painterTc.setPen( QPen( QColor( 70, 70, 70 ) ) );
    painterTc.setFont( font );
    painterTc.drawText( 70, 14, 380, 40, 0, QString( "88 : 88 : 88 . 88" ) );

    //Foreground
    painterTc.setPen( QPen( QColor( 220, 220, 220 ) ) );
    painterTc.setFont( font );
    painterTc.drawText( 70, 14, 380, 40, 0, QString( "%1 : %2 : %3 . %4" )
                      .arg( hours, 2, 10, QChar('0') )
                      .arg( minutes, 2, 10, QChar('0') )
                      .arg( seconds, 2, 10, QChar('0') )
                      .arg( frame, 2, 10, QChar('0') ) );

    return *m_tcImage;
}